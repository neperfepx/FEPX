! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE MATRIX_OPERATIONS_MOD
!
! This module contains a variety of matrix operations and other utility
!   functions.
!
! Contains subroutines:
! CALC_ELVOL: Calculates element volumes.
! DETERMINANT_GRN: Calculates determinants of an array of 3x3 matrices
! FIND_INDICES: Find the members of vector equal to target, return val and ind.
! GEN_MATRIX_VECTOR_MULT: Matrix times vector
! GEN_MATRIX_MULT: Matrix times matrix
! INVERT5X5: Inverts an array of 5x5 matrices.
! LATTICE_DEFORM: Change coordinates from sample to lattice (for 5-vectors).
! LATTICE_DEFORM_SER: Serial version (one element).
! LATTICE_SPIN: Convert a skew 3-vector from sample to lattice coordinates.
! LATTICE_SPIN_SER: Serial version (one element).
! MAT_VEC_SKEW: Convert array of skew matrices to array of 3-vectors.
! MAT_VEC_SKEW_SER: Serial version (one element).
! MAT_VEC_SYMM: Convert array of 3x3 symmetric matrices to array of 5-vectors.
! MAT_VEC_SYMM_SER: Serial version (one element).
! MAT_X_MAT3: Matrix multiplication for arrays of 3x3 matrices. (c = a*b)
! MAT_X_MAT5: Matrix multiplication for arrays of 5x5 matrices.
! MAT_X_MATT3: Matrix multiplication by transpose. (3x3)
! MAT_X_MATT5: Matrix multiplication by transpose for arrays of 5x5 matrices.
! MAT_X_MATS3: Multiply array of 3x3 matrices by a fixed 3x3 matrix, tranposed.
! MAT_X_VEC5: Multiply array of matrices times array of vectors. (5 dim)
! MAT_X_VEC5_SER: Serial version (one element).
! MATRIX_VEC_MULT: Matrix-vector multiplication routine (RC, 2017)
! MATT_X_MAT3: Matrix multiplication by transpose (first matrix T). (3x3)
! NORM_VEC: Compute 2-norm for array of 5-vectors.
! ROT_MAT_SKEW: Construct 3x3 rotation matrix acting on skew-matrix 3-vectors.
! ROT_MAT_SKEW_SER: Serial version (one element).
! ROT_MAT_SYMM: Construct 5x5 rotation matrix acting on 5-vectors.
! ROT_MAT_SYMM_SER: Serial version (one element).
! SOLVE_LIN_SYS_3: Solve a linear system of three equations. Used for triaxial
!   loading.
! SPARSE_MATVEC_EBE: Matrix times vector, element by element
! SYMM_VGR: Compute symmetric part of an array of velocity gradients.
! SYMM_VGR_SER: Serial version (one element).
! SKEW_VGR: Compute skew part of an array of velocity gradients.
! TENSOR3DCOMPOSE: Form matrix from deviatoric, skew, or spherical parts. If no
!   parts are passed, the matix is zeroed.
! TESNRO3DDECOMPOSE: Decompose matrix into deviatoric and spherical parts.
! VEC5_VEC6: Convert 5-vector to 6-vector of symmetric matrix.
! VEC_D_VEC5: Multiply diagonal matrix times array of vectors. (5 dim)
! VEC_MAT_SKEW: Convert array of 3-vectors to array of skew matrices.
! VEC_MAT_SKEW_GRN: Also considers legacy grains per element.
! VEC_MAT_SYMM: Convert 5-vector to symmetric matrix.
! VEC_MAT_SYMM_GRN: Also considers legacy grains per element.
! VEC_MAT_SYMM_SER: Serial version (one element).
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND
!
! From libfepx:
!
USE DIMENSIONS_MOD
USE QUADRATURE_MOD
USE READ_INPUT_MOD
USE SHAPE_3D_MOD
!
! From libparallel:
!
USE GATHER_SCATTER_MOD
!
IMPLICIT NONE
!
! Parameters (some private, some public)
!
INTEGER, PARAMETER:: DECOMP_MPSIM = 0
INTEGER, PARAMETER:: DECOMP_FEMEVPS = 1
INTEGER, PRIVATE :: DECOMP_DFLT = DECOMP_MPSIM
!
! Constants (all private)
!
REAL(RK), PARAMETER, PRIVATE :: SQ2_I = 1.0D0 / DSQRT(2.0D0)
REAL(RK), PARAMETER, PRIVATE :: SQ6_I = 1.0D0 / DSQRT(6.0D0)
REAL(RK), PARAMETER, PRIVATE :: TWOSQ6_I = 2.0D0 * DSQRT(6.0D0)
!
CONTAINS 
    !
    SUBROUTINE CALC_ELVOL(ELVOL, ECOORDS)
    !
    ! Calculates element volumes
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! ELVOL: Array of elemental volumes
    ! ECOORDS: Coordinates of elemental nodal points
    !
    REAL(RK), INTENT(OUT) :: ELVOL(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: ECOORDS(0:KDIM1, EL_SUB1:EL_SUP1)
    !
    ! Locals:
    ! I: Looping index
    ! All others: Intermediary calculations
    !
    REAL(RK) :: LOC0
    REAL(RK) :: LOC1
    REAL(RK) :: LOC2
    REAL(RK) :: WT(EL_SUB1:EL_SUP1)
    REAL(RK) :: DET(EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDX(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDY(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDZ(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK) :: S11(EL_SUB1:EL_SUP1)
    REAL(RK) :: S12(EL_SUB1:EL_SUP1)
    REAL(RK) :: S13(EL_SUB1:EL_SUP1)
    REAL(RK) :: S21(EL_SUB1:EL_SUP1)
    REAL(RK) :: S22(EL_SUB1:EL_SUP1)
    REAL(RK) :: S23(EL_SUB1:EL_SUP1)
    REAL(RK) :: S31(EL_SUB1:EL_SUP1)
    REAL(RK) :: S32(EL_SUB1:EL_SUP1)
    REAL(RK) :: S33(EL_SUB1:EL_SUP1)
    INTEGER ::  I
    !
    !---------------------------------------------------------------------------
    !
    ELVOL = 0.0D0
    !
    DO I = 0, NQPT1
        !
        LOC0 = QPLOC(0, I)
        LOC1 = QPLOC(1, I)
        LOC2 = QPLOC(2, I)
        WT = WTQP(0, I)
        !
        CALL SFDER_HPAR(LOC0, LOC1, LOC2, ECOORDS, DNDX, DNDY, DNDZ, DET, S11, &
            & S12, S13, S21, S22, S23, S31, S32, S33)
        !
        ELVOL = ELVOL + DET * WT
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE CALC_ELVOL
    !
    !===========================================================================
    !
    SUBROUTINE DETERMINANT_GRN(TENSOR, DETERM, N, M)
    !
    ! Calculates determinants of an array of 3x3 matrices.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! TENSOR: Array of 3x3 matrices
    ! DETERM: Array of scalars (determinant of tensor)
    ! N: Number of grains/element (legacy, always 0:0)
    ! M: Number of elements
    !
    REAL(RK), INTENT(IN) :: TENSOR(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: DETERM(0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    ! I: Looping index
    ! DETxx: Intermediary calculations
    !
    INTEGER :: I
    REAL(RK) :: DET11(0:(N - 1), 0:(M - 1))
    REAL(RK) :: DET12(0:(N - 1), 0:(M - 1))
    REAL(RK) :: DET13(0:(N - 1), 0:(M - 1))
    REAL(RK) :: DET21(0:(N - 1), 0:(M - 1))
    REAL(RK) :: DET22(0:(N - 1), 0:(M - 1))
    REAL(RK) :: DET23(0:(N - 1), 0:(M - 1))
    !
    !---------------------------------------------------------------------------
    !
    DET11 = TENSOR(0, 0, :, :)* TENSOR(1, 1, :, :)* TENSOR(2, 2, :, :)
    DET12 = TENSOR(0, 1, :, :)* TENSOR(1, 2, :, :)* TENSOR(2, 0, :, :)
    DET13 = TENSOR(0, 2, :, :)* TENSOR(1, 0, :, :)* TENSOR(2, 1, :, :)
    DET21 = TENSOR(0, 0, :, :)* TENSOR(1, 2, :, :)* TENSOR(2, 1, :, :)
    DET22 = TENSOR(0, 1, :, :)* TENSOR(1, 0, :, :)* TENSOR(2, 2, :, :)
    DET23 = TENSOR(0, 2, :, :)* TENSOR(1, 1, :, :)* TENSOR(2, 0, :, :)
    !
    DETERM = DET11 + DET12 + DET13 - DET21 - DET22 - DET23
    !
    RETURN
    !
    END SUBROUTINE DETERMINANT_GRN
    !
    !===========================================================================
    !
    SUBROUTINE FIND_INDICES(NUMIND, TARGET, VECTOR, INDICES)
    !
    ! Find the members of vector (integer) which are equal to target (integer),
    !   and return indices and numind
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    !
    INTEGER :: NUMIND
    INTEGER :: TARGET
    INTEGER :: VECTOR(:)
    INTEGER, POINTER :: INDICES(:)
    !
    ! Locals:
    !
    INTEGER :: K
    INTEGER :: LB
    INTEGER :: UB
    INTEGER :: I
    !
    !---------------------------------------------------------------------------
    !
    ! Tito Marin (legacy):
    !
    NULLIFY(INDICES)
    !
    ! See if indices has been allocated. If so, deallocate.
    !
    IF (ASSOCIATED(INDICES)) THEN
        !
        DEALLOCATE(INDICES)
        !
    END IF
    !
    ! Find how many indices there are.
    !
    NUMIND = COUNT(VECTOR .EQ. TARGET)
    !
    ! Allocate indices to be numind long.
    !
    ALLOCATE(INDICES(1:NUMIND))
    !
    ! Find which indices of vector equal target, and save them in indices.
    !
    K = 0
    LB = LBOUND(VECTOR, 1) ! Always equals 1
    UB = UBOUND(VECTOR, 1)
    !
    DO I = LB, UB
        !
        IF (VECTOR(I) .EQ. TARGET) THEN
            !
            K = K + 1
            INDICES(K) = I
            !
        END IF
        !
    END DO
    !
    ! 0-based
    !
    INDICES = INDICES - 1
    !
    RETURN
    !
    END SUBROUTINE FIND_INDICES
    !
    !===========================================================================
    !
    SUBROUTINE GEN_MATRIX_VECTOR_MULT(Y, A, X, I1, I2, I3, I4, IER)
    !
    !  Array of matrices times array of vectors: y(i) = A(i)*x(i)
    !
    ! Note that the arrays x,y,a are all assumed shape arrays and that they are
    !   therefore dimensioned 1:p, 1:q, etc., even though they may be
    !   dimensioned differently in the calling routine.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! Y:
    ! A:
    ! X:
    ! I1:
    ! I2:
    ! I3:
    ! I4:
    ! IER:
    !
    REAL(RK) :: Y(:,:)
    REAL(RK) :: A(:,:,:)
    REAL(RK) :: X(:,:)
    INTEGER :: I1
    INTEGER :: I2
    INTEGER :: I3
    INTEGER :: I4
    INTEGER :: IER
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: K
    INTEGER :: AD1
    INTEGER :: AD2
    INTEGER :: AD3
    !
    !---------------------------------------------------------------------------
    !
    AD1 = UBOUND(A, 1)
    AD2 = UBOUND(A, 2)
    AD3 = UBOUND(A, 3)
    !
    DO K = 1, AD3
        !
        DO J = 1, AD2
            !
            !DO I = 1, AD1
            Y(:, K) = Y(:, K) + A(:, J, K) * X(J, K)
            !END DO
            !
        END DO
        !
    END DO
    !
    IER = 0
    !
    RETURN
    !
    END SUBROUTINE GEN_MATRIX_VECTOR_MULT
    !
    !===========================================================================
    !
    SUBROUTINE GEN_MATRIX_MULT(C, A, B, I1, I2, IER)
    !
    ! Accumulate Array of matrices times array of matrices: c(i) = a(i) * b(i)
    !
    ! Note: c is not zeroed in this subroutine, so the values are accumulated;
    !   in all the calling routines, however, the c array is zeroed before
    !   calling this routine.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! C: Array of result matrices (c=a*b)
    ! A, B: Array of input matrices
    ! I1, I2: Not used (MPK: So why are they here?)
    ! IER: Return status (0 = okay)
    !
    REAL(RK) :: C(:,:,:)
    REAL(RK) :: A(:,:,:)
    REAL(RK) :: B(:,:,:)
    INTEGER :: I1
    INTEGER :: I2
    INTEGER :: IER
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: K
    INTEGER :: M
    INTEGER :: LDA
    INTEGER :: LDB
    INTEGER :: LDC
    INTEGER :: LTA
    INTEGER :: LTB
    INTEGER :: N
    !
    !---------------------------------------------------------------------------
    !
    LDA = UBOUND(A, 1)
    LDB = UBOUND(B, 1)
    LDC = UBOUND(C, 1)
    LTA = UBOUND(A, 2)
    LTB = UBOUND(B, 2)
    N = UBOUND(A, 3)
    !
    DO M = 1, N
        !
        DO J = 1, LTB
            !
            DO K = 1, LTA
                !
!               DO I = 1,LDC
                C(:, J, M) = C(:, J, M) + A(:, K, M) * B(K, J, M)
!               END DO
                !
            END DO
            !
        END DO
        !
    END DO
    !
    IER = 0
    !
    RETURN
    !
    END SUBROUTINE GEN_MATRIX_MULT
    !
    !===========================================================================
    !
    SUBROUTINE INVERT5X5(A, N, M)
    !
    ! Invert array of 5x5 matrices.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    ! A: Array of matrices on input and array of inverses on output
    ! N: Number of grains/element (legacy, always 0:0)
    ! M: Number of elements
    !
    REAL(RK), INTENT(INOUT) :: A(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    ! I, J: Looping indices
    ! Axx, AIxx: Intermediary calculations
    !
    INTEGER :: I
    INTEGER :: J
    REAL(RK) :: A11(0:(N - 1), 0:(M - 1))
    REAL(RK) :: A21(0:(N - 1), 0:(M - 1))
    REAL(RK) :: A22(0:(N - 1), 0:(M - 1))
    REAL(RK) :: A31(0:(N - 1), 0:(M - 1))
    REAL(RK) :: A32(0:(N - 1), 0:(M - 1))
    REAL(RK) :: A33(0:(N - 1), 0:(M - 1))
    REAL(RK) :: A41(0:(N - 1), 0:(M - 1))
    REAL(RK) :: A42(0:(N - 1), 0:(M - 1))
    REAL(RK) :: A43(0:(N - 1), 0:(M - 1))
    REAL(RK) :: A44(0:(N - 1), 0:(M - 1))
    REAL(RK) :: A51(0:(N - 1), 0:(M - 1))
    REAL(RK) :: A52(0:(N - 1), 0:(M - 1))
    REAL(RK) :: A53(0:(N - 1), 0:(M - 1))
    REAL(RK) :: A54(0:(N - 1), 0:(M - 1))
    REAL(RK) :: A55(0:(N - 1), 0:(M - 1))
    REAL(RK) :: AI11(0:(N - 1), 0:(M - 1))
    REAL(RK) :: AI21(0:(N - 1), 0:(M - 1))
    REAL(RK) :: AI22(0:(N - 1), 0:(M - 1))
    REAL(RK) :: AI31(0:(N - 1), 0:(M - 1))
    REAL(RK) :: AI32(0:(N - 1), 0:(M - 1))
    REAL(RK) :: AI33(0:(N - 1), 0:(M - 1))
    REAL(RK) :: AI41(0:(N - 1), 0:(M - 1))
    REAL(RK) :: AI42(0:(N - 1), 0:(M - 1))
    REAL(RK) :: AI43(0:(N - 1), 0:(M - 1))
    REAL(RK) :: AI44(0:(N - 1), 0:(M - 1))
    REAL(RK) :: AI51(0:(N - 1), 0:(M - 1))
    REAL(RK) :: AI52(0:(N - 1), 0:(M - 1))
    REAL(RK) :: AI53(0:(N - 1), 0:(M - 1))
    REAL(RK) :: AI54(0:(N - 1), 0:(M - 1))
    REAL(RK) :: AI55(0:(N - 1), 0:(M - 1))
    !
    !---------------------------------------------------------------------------
    !
    A11 = A(0, 0, :, :)
    A21 = A(1, 0, :, :)
    A22 = A(1, 1, :, :)
    A31 = A(2, 0, :, :)
    A32 = A(2, 1, :, :)
    A33 = A(2, 2, :, :)
    A41 = A(3, 0, :, :)
    A42 = A(3, 1, :, :)
    A43 = A(3, 2, :, :)
    A44 = A(3, 3, :, :)
    A51 = A(4, 0, :, :)
    A52 = A(4, 1, :, :)
    A53 = A(4, 2, :, :)
    A54 = A(4, 3, :, :)
    A55 = A(4, 4, :, :)
    !
    ! ** A = LDL'.
    ! j = 1
    !
    AI55 = 1.0D0 / A11
    A21 = A21 * AI55
    A31 = A31 * AI55
    A41 = A41 * AI55
    A51 = A51 * AI55
    !
    ! j = 2
    !
    AI11 = A21 * A11
    A22 = A22 - A21 * AI11
    !
    A32 = A32 - A31 * AI11
    A42 = A42 - A41 * AI11
    A52 = A52 - A51 * AI11
    AI55 = 1.0D0 / A22
    A32 = A32 * AI55
    A42 = A42 * AI55
    A52 = A52 * AI55
    !
    ! j = 3
    !
    AI11 = A31 * A11
    AI22 = A32 * A22
    AI55 = A31 * AI11 + A32 * AI22
    !
    A33 = A33 - AI55
    !
    AI55 = 1.0D0 / A33
    A43 = A43 - A41 * AI11 - A42 * AI22
    A53 = A53 - A51 * AI11 - A52 * AI22
    A43 = A43 * AI55
    A53 = A53 * AI55
    !
    ! j = 4
    !
    AI11 = A41 * A11
    AI22 = A42 * A22
    AI33 = A43 * A33
    AI55 = A41 * AI11 + A42 * AI22 + A43 * AI33
    !
    A44 = A44 - AI55
    !
    A54 = A54 - A51 * AI11 - A52 * AI22 - A53 * AI33
    A54 = A54 / A44
    !
    ! j = 5
    !
    AI11 = A51 * A11
    AI22 = A52 * A22
    AI33 = A53 * A33
    AI44 = A54 * A44
    AI55 = A51 * AI11 + A52 * AI22 + A53 * AI33 + A54 * AI44
    !
    A55 = A55 - AI55
    !
    ! Column 1 of inverse
    ! Ly = b
    !
    AI21 = - A21
    AI31 = - A31 - A32 * AI21
    AI41 = - A41 - A42 * AI21 - A43 * AI31
    AI51 = - A51 - A52 * AI21 - A53 * AI31 - A54 * AI41
    !
    ! Dz = y
    !
    AI11 = 1.0D0 / A11
    AI21 = AI21 / A22
    AI31 = AI31 / A33
    AI41 = AI41 / A44
    AI51 = AI51 / A55
    !
    ! L'x = z
    !
    AI41 = AI41 - A54 * AI51
    AI31 = AI31 - A43 * AI41 - A53 * AI51
    AI21 = AI21 - A32 * AI31 - A42 * AI41 - A52 * AI51
    AI11 = AI11 - A21 * AI21 - A31 * AI31 - A41 * AI41 - A51 * AI51
    !
    ! Column 2 of inverse
    ! Ly = b
    !
    AI32 = - A32
    AI42 = - A42 - A43 * AI32
    AI52 = - A52 - A53 * AI32 - A54 * AI42
    !
    ! Dz = y
    !
    AI22 = 1.0D0 / A22
    AI32 = AI32 / A33
    AI42 = AI42 / A44
    AI52 = AI52 / A55
    !
    ! L'x = z
    !
    AI42 = AI42 - A54 * AI52
    AI32 = AI32 - A43 * AI42 - A53 * AI52
    AI22 = AI22 - A32 * AI32 - A42 * AI42 - A52 * AI52
    !
    ! Column 3 of inverse
    ! Ly = b
    !
    AI43 = - A43
    AI53 = - A53 - A54 * AI43
    !
    ! Dz = y
    !
    AI33 = 1.0D0 / A33
    AI43 = AI43 / A44
    AI53 = AI53 / A55
    !
    ! L'x = z
    !
    AI43 = AI43 - A54 * AI53
    AI33 = AI33 - A43 * AI43 - A53 * AI53
    !
    ! Column 4 of inverse
    ! Ly = b
    !
    AI54 = - A54
    !
    ! Dz = y
    !
    AI44 = 1.0D0 / A44
    AI54 = AI54 / A55
    !
    ! L'x = z
    !
    AI44 = AI44 - A54 * AI54
    !
    ! Column 5 of inverse
    ! Dz = y
    !
    AI55 = 1.0D0 / A55
    !
    ! Recover Array
    !
    A(0, 0, :, :) = AI11
    A(1, 0, :, :) = AI21
    A(2, 0, :, :) = AI31
    A(3, 0, :, :) = AI41
    A(4, 0, :, :) = AI51
    A(1, 1, :, :) = AI22
    A(2, 1, :, :) = AI32
    A(3, 1, :, :) = AI42
    A(4, 1, :, :) = AI52
    A(2, 2, :, :) = AI33
    A(3, 2, :, :) = AI43
    A(4, 2, :, :) = AI53
    A(3, 3, :, :) = AI44
    A(4, 3, :, :) = AI54
    A(4, 4, :, :) = AI55
    !
    DO I = 0, TVEC1
        !
        DO J = (I + 1), TVEC1
            !
            A(I, J, :, :) = A(J, I, :, :)
            !
        END DO
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE INVERT5X5
    !
    !===========================================================================
    !
    SUBROUTINE LATTICE_DEFORM(QR5X5, VEC, VEC_LAT, N, M)
    !
    ! Change coordinates from sample to lattice (for 5-vectors).
    !
    ! Rotate{VEC}_SAM to {VEC}_LAT via:
    ! {VEC_LAT} = [QR5X5]' * {VEC}
    !
    !----------------------------------------------------------------------
    !
    ! Arguments:
    ! QR5X5: Array of rotation matrices from lattice to sample
    ! VEC: Array of 5-vectors in sample basis
    ! VEC_LAT: Array of 5-vectors in lattice basis
    ! N: Number of grains/element (legacy, always 0:0)
    ! M: Number of elements
    !
    REAL(RK), INTENT(IN) :: QR5X5(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: VEC(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: VEC_LAT(0:TVEC1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    ! I, J: Looping indices
    !
    INTEGER :: I
    INTEGER :: J
    !
    !----------------------------------------------------------------------
    !
    ! Rotate {vec}_sm -> {vec}_lat : {v_lat} = [QR5X5]' {v_sm}
    !
    VEC_LAT = 0.0D0
    !
    DO I = 0, TVEC1
        !
        DO J = 0, TVEC1
            !
            VEC_LAT(I, :, :) = VEC_LAT(I, :, :) + QR5X5(J, I, :, :) * &
                & VEC(J, :, :)
            !
        END DO
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE LATTICE_DEFORM
    !
    !===========================================================================
    !
    SUBROUTINE LATTICE_DEFORM_SER(QR5X5, VEC, VEC_LAT)
    !
    ! Change coordinates from sample to lattice (for 5-vectors).
    !
    ! Rotate {VEC}_SAM to {VEC}_LAT via:
    ! {VEC_LAT} = [QR5X5]' * {VEC}
    !
    !----------------------------------------------------------------------
    !
    ! Arguments:
    ! QR5X5: Rotation matrix from lattice to sample
    ! VEC: 5-vector in sample basis
    ! VEC_LAT: 5-vector in lattice basis
    !
    REAL(RK), INTENT(IN) :: QR5X5(0:TVEC1, 0:TVEC1)
    REAL(RK), INTENT(IN) :: VEC(0:TVEC1)
    REAL(RK), INTENT(OUT) :: VEC_LAT(0:TVEC1)
    !
    ! Locals:
    ! I, J: Looping indices
    !
    INTEGER :: I
    INTEGER :: J
    !
    !----------------------------------------------------------------------
    !
    VEC_LAT = 0.0D0
    !
    DO I = 0, TVEC1
        !
        DO J = 0, TVEC1
            !
            VEC_LAT(I) = VEC_LAT(I) + QR5X5(J, I) * VEC(J)
            !
        END DO
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE LATTICE_DEFORM_SER
    !
    !===========================================================================
    !
    SUBROUTINE LATTICE_SPIN(QR3X3, W_VEC, W_VEC_LAT, N, M)
    !
    ! Convert a skew 3-vector from sample to lattice coordinates.
    !
    ! Rotate {W_VEC}_SAM to {W_VEC}_LAT via:
    ! {W_VEC_LAT} = [QR3X3]' * {W_VEC}
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! QR3X3: Array of rotation matrices from lattice to sample
    ! W_VEC: Array of skew 3-vectors in sample basis
    ! W_VEC_LAT: Array of skew 3-vectors in lattice basis
    ! N: Number of grains/element (legacy, always 0:0)
    ! M: Number of elements
    !
    REAL(RK), INTENT(IN) :: QR3X3(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: W_VEC(0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: W_VEC_LAT(0:DIMS1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    ! I, J: Looping indices
    !
    INTEGER :: I
    INTEGER :: J
    !
    !---------------------------------------------------------------------------
    !
    W_VEC_LAT = 0.0D0
    !
    DO I = 0, DIMS1
        !
        DO J = 0, DIMS1
            !
            W_VEC_LAT(I, :, :) = W_VEC_LAT(I, :, :) + QR3X3(J, I, :, :) * &
                & W_VEC(J, :, :)
            !
        END DO
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE LATTICE_SPIN
    !
    !===========================================================================
    !
    SUBROUTINE LATTICE_SPIN_SER(QR3X3, W_VEC, W_VEC_LAT)
    !
    ! Convert a skew 3-vector from sample to lattice coordinates.
    !
    ! Rotate {W_VEC}_SAM to {W_VEC}_LAT via:
    ! {W_VEC_LAT} = [QR3X3]' * {W_VEC}
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! QR3X3: Rotation matrix from lattice to sample (to be transposed)
    ! W_VEC: Skew 3-vector in sample basis
    ! W_VEC_LAT: Skew 3-vector in lattice basis
    !
    REAL(RK), INTENT(IN) :: QR3X3(0:DIMS1, 0:DIMS1)
    REAL(RK), INTENT(IN) :: W_VEC(0:DIMS1)
    REAL(RK), INTENT(OUT) :: W_VEC_LAT(0:DIMS1)
    !
    ! Locals:
    ! I, J: Looping indices
    !
    INTEGER :: I
    INTEGER :: J
    !
    !---------------------------------------------------------------------------
    !
    W_VEC_LAT = 0.0D0
    !
    DO I = 0, DIMS1
        !
        DO J = 0, DIMS1
            !
            W_VEC_LAT(I) = W_VEC_LAT(I) + QR3X3(J, I) * W_VEC(J)
            !
        END DO
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE LATTICE_SPIN_SER
    !
    !===========================================================================
    !
    SUBROUTINE MAT_VEC_SKEW(MAT, VEC, M)
    !
    ! Convert array of skew matrices to array of 3-vectors.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! MAT: Array of skew matrices
    ! VEC: Array of 3-vectors
    ! M: Number of elements
    !
    REAL(RK), INTENT(IN) :: MAT(0:DIMS1, 0:DIMS1, 0:(M - 1))
    REAL(RK), INTENT(OUT) :: VEC(0:DIMS1, 0:(M - 1))
    INTEGER, INTENT(IN) :: M
    !
    !---------------------------------------------------------------------------
    !
    VEC = 0.0D0
    !
    VEC(0, :) = MAT(1, 0, :)
    VEC(1, :) = MAT(2, 0, :)
    VEC(2, :) = MAT(2, 1, :)
    !
    RETURN
    !
    END SUBROUTINE MAT_VEC_SKEW
    !
    !===========================================================================
    !
    SUBROUTINE MAT_VEC_SKEW_SER(MAT, VEC)
    !
    ! Convert skew 3x3 matrix to 3-vector.
    !
    !----------------------------------------------------------------------
    !
    ! Arguments:
    ! MAT: Skew matrix
    ! VEC: 3-vector
    !
    REAL(RK), INTENT(IN) :: MAT(0:DIMS1, 0:DIMS1)
    REAL(RK), INTENT(OUT) :: VEC(0:DIMS1)
    !
    !----------------------------------------------------------------------
    !
    VEC = 0.0D0
    !
    VEC(0) = MAT(1, 0)
    VEC(1) = MAT(2, 0)
    VEC(2) = MAT(2, 1)
    !
    RETURN
    !
    END SUBROUTINE MAT_VEC_SKEW_SER
    !
    !===========================================================================
    !
    SUBROUTINE MAT_VEC_SYMM(MAT, VEC, M)
    !
    ! Convert array of 3x3 symmetric matrices to array of 5-vectors.
    !
    !----------------------------------------------------------------------
    !
    ! Arguments:
    ! MAT: Array of symmetric matrices
    ! VEC: Array of 5-vectors
    ! M: Number of elements
    !
    REAL(RK), INTENT(IN) :: MAT(0:DIMS1, 0:DIMS1, 0:(M - 1))
    REAL(RK), INTENT(OUT) :: VEC(0:TVEC1, 0:(M - 1))
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    ! SQR2, SQR32: Local parameters
    !
    REAL(RK) :: SQR2
    REAL(RK) :: SQR32
    !
    !----------------------------------------------------------------------
    !
    SQR2 = DSQRT(2.0D0)
    SQR32 = DSQRT(1.5D0)
    VEC = 0.0D0
    !
    VEC(0, :) = (MAT(0, 0, :) - MAT(1, 1, :)) / SQR2
    VEC(1, :) = MAT(2, 2, :) * SQR32
    VEC(2, :) = MAT(1, 0, :) * SQR2
    VEC(3, :) = MAT(2, 0, :) * SQR2
    VEC(4, :) = MAT(2, 1, :) * SQR2
    !
    RETURN
    !
    END SUBROUTINE MAT_VEC_SYMM
    !
    !===========================================================================
    !
    SUBROUTINE MAT_VEC_SYMM_SER(MAT, VEC)
    !
    ! Convert 3x3 symmetric matrix to 5-vector.
    !
    !----------------------------------------------------------------------
    !
    ! Arguments:
    ! MAT: Symmetric matrix
    ! VEC: 5-vector
    !
    REAL(RK), INTENT(IN) :: MAT(0:DIMS1, 0:DIMS1)
    REAL(RK), INTENT(OUT) :: VEC(0:TVEC1)
    !
    ! Locals:
    ! SQR2, SQR32: Local parameters
    !
    REAL(RK) :: SQR2
    REAL(RK) :: SQR32
    !
    !----------------------------------------------------------------------
    !
    SQR2 = DSQRT(2.0D0)
    SQR32 = DSQRT(1.5D0)
    VEC = 0.0D0
    !
    VEC(0) = (MAT(0, 0) - MAT(1, 1)) / SQR2
    VEC(1) = MAT(2, 2) * SQR32
    VEC(2) = MAT(1, 0) * SQR2
    VEC(3) = MAT(2, 0) * SQR2
    VEC(4) = MAT(2, 1) * SQR2
    !
    RETURN
    !
    END SUBROUTINE MAT_VEC_SYMM_SER
    !
    !===========================================================================
    !
    SUBROUTINE MAT_X_MAT3(A, B, C, N, M)
    !
    ! Matrix multiplication for arrays of 3x3 matrices. (c = a*b)
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! A, B: Arrays of input matrices
    ! C: Array of output matrices (a*b)
    ! N: Number of grains/element (legacy, always 0:0)
    ! M: Number of elements
    !
    REAL(RK), INTENT(IN) :: A(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: B(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: C(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    ! I, J, K: Looping indices
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: K
    !
    !---------------------------------------------------------------------------
    !
    C = 0.0D0
    !
    DO J = 0, DIMS1
        !
        DO K = 0, DIMS1
            !
            DO I = 0, DIMS1
                !
                C(I, J, :, :) = C(I, J, :, :) + A(I, K, :, :) * B(K, J, :, :)
                !
            END DO
            !
        END DO
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE MAT_X_MAT3
    !
    !===========================================================================
    !
    SUBROUTINE MAT_X_MAT5(A, B, C, N, M)
    !
    ! Matrix multiplication for arrays of 5x5 matrices.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! A, B: Arrays of 5x5 matrices (input)
    ! C: Array of 5x5 matrices (c=a*b) (output)
    ! N: Number of grains/element (legacy, always 0:0)
    ! M: Number of elements
    !
    REAL(RK), INTENT(IN) :: A(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: B(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: C(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    ! I, J, K: Looping indices
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: K
    !
    !----------------------------------------------------------------------
    !
    C = 0.0D0
    !
    DO J = 0, TVEC1
        !
        DO K = 0, TVEC1
            !
            DO I = 0, TVEC1
                !
                C(I, J, :, :) = C(I, J, :, :) + A(I, K, :, :) * B(K, J, :, :)
                !
            END DO
            !
        END DO
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE MAT_X_MAT5
    !
    !===========================================================================
    !
    SUBROUTINE MAT_X_MATT3(A, B, C, N, M)
    !
    ! Matrix multiplication by transpose. (3x3)
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! A, B: Arrays of 3x3 matrices
    ! C: Array of 3x3 matrices, c=a*b^t
    ! N: Number of grains/element (legacy, always 0:0)
    ! M: Number of elements
    !
    REAL(RK), INTENT(IN) :: A(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: B(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: C(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    ! I, J, K: Looping indices
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: K
    !
    !---------------------------------------------------------------------------
    !
    C = 0.0D0
    !
    DO J = 0, DIMS1
        !
        DO K = 0, DIMS1
            !
            DO I = 0, DIMS1
                !
                C(I, J, :, :) = C(I, J, :, :) + A(I, K, :, :) * B(J, K, :, :)
                !
            END DO
            !
        END DO
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE MAT_X_MATT3
    !
    !===========================================================================
    !
    SUBROUTINE MAT_X_MATT5(A, B, C, N, M)
    !
    ! Matrix multiplication by transpose for arrays of 5x5 matrices.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! A, B: Arrays of 5x5 matrices (input)
    ! C: Array of 5x5 matrices; c=a*b^t (output)
    ! N: Number of grains/element (legacy, always 0:0)
    ! M: Number of elements
    !
    REAL(RK), INTENT(IN) :: A(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: B(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: C(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN):: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    ! I, J, K: Looping indices
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: K
    !
    !---------------------------------------------------------------------------
    !
    C = 0.0D0
    !
    DO K = 0, TVEC1
        !
        DO J = 0, TVEC1
            !
            DO I = 0, TVEC1
                !
                C(I, J, :, :) = C(I, J, :, :) + A(I, K, :, :) * B(J, K, :, :)
                !
            END DO
            !
        END DO
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE MAT_X_MATT5
    !
    !===========================================================================
    !
    SUBROUTINE MAT_X_MATS3(A, B, C, N, M)
    !
    !     Multiply array of 3x3 matrices by a fixed 3x3 matrix, tranposed.
    !
    !----------------------------------------------------------------------
    !
    ! Arguments:
    ! A: Array of matrices to be multiplied
    ! B: Fixed 3x3 matrix
    ! C: Array of 3x3 matrices, c=a*b^t
    ! N: Number of grains/element (legacy, always 0:0)
    ! M: Number of elements
    !
    REAL(RK), INTENT(IN) :: A(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: B(0:DIMS1, 0:DIMS1)
    REAL(RK), INTENT(OUT) :: C(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    ! I, J, K: Looping indices
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: K
    !
    !----------------------------------------------------------------------
    !
    C = 0.0D0
    !
    DO K = 0, DIMS1
        !
        DO J = 0, DIMS1
            !
            DO I = 0, DIMS1
                !
                C(I, J, :, :) = C(I, J, :, :) + A(I, K, :, :) * B(J, K)
                !
            END DO
            !
        END DO
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE MAT_X_MATS3
    !
    !===========================================================================
    !
    SUBROUTINE MAT_X_VEC5(MATRIX, VECTOR, PRODUCT, N, M)
    !
    ! Multiply array of matrices times array of vectors. (5 dim)
    !
    !----------------------------------------------------------------------
    !
    ! Arguments:
    ! MATRIX: Array of matrices (5x5)
    ! VECTOR: Array of 5-vectors
    ! PRODUCT: Array of 5-vectors; product=matrix * vector
    ! N: Number of grains/element (legacy, always 0:0)
    ! M: Number of elements
    !
    REAL(RK), INTENT(IN) :: MATRIX(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: VECTOR(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: PRODUCT(0:TVEC1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    ! I, J: Looping indices
    !
    INTEGER :: I
    INTEGER :: J
    !
    !----------------------------------------------------------------------
    !
    PRODUCT = 0.0D0
    !
    DO J = 0, TVEC1
        !
        DO I = 0, TVEC1
            !
            PRODUCT(I, :, :) = PRODUCT(I, :, :) + MATRIX(I, J, :, :) * &
                & VECTOR(J, :, :)
            !
        END DO
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE MAT_X_VEC5
    !
    !===========================================================================
    !
    SUBROUTINE MAT_X_VEC5_SER(MATRIX, VECTOR, PRODUCT)
    !
    ! Multiply array of matrices times array of vectors. (5 dim)
    !
    !----------------------------------------------------------------------
    !
    ! Arguments:
    ! MATRIX: Matrix (5x5)
    ! VECTOR: 5-vectors
    ! PRODUCT: 5-vector; product=matrix * vector
    !
    REAL(RK), INTENT(IN) :: MATRIX(0:TVEC1, 0:TVEC1)
    REAL(RK), INTENT(IN) :: VECTOR(0:TVEC1)
    REAL(RK), INTENT(OUT) :: PRODUCT(0:TVEC1)
    !
    ! Locals:
    ! I, J: Looping indices
    !
    INTEGER :: I
    INTEGER :: J
    !
    !----------------------------------------------------------------------
    !
    PRODUCT = 0.0D0
    !
    DO J = 0, TVEC1
        !
        DO I = 0, TVEC1
            !
            PRODUCT(I) = PRODUCT(I) + MATRIX(I, J) * VECTOR(J)
            !
        END DO
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE MAT_X_VEC5_SER
    !
    !===========================================================================
    !
    SUBROUTINE MATRIX_VEC_MULT(M, V, MV, DIM)
    !
    ! A simple matrix vector multiplication routine that takes advantage of the
    !   compiler's vectorization optimizations. It ends up being faster than
    !   using the built in Fortran matmul method (RC 2017).
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! M: Array of matrices
    ! V: Array of vectors
    ! MV: Array of products of M*V
    ! DIM:
    !
    REAL(RK), INTENT(IN) :: M(0:DIM-1,0:DIM-1)
    REAL(RK), INTENT(IN) :: V(0:DIM-1)
    REAL(RK), INTENT(OUT) :: MV(0:DIM-1)
    INTEGER, INTENT(IN) :: DIM
    !
    ! Locals:
    !
    INTEGER :: I
    !
    !---------------------------------------------------------------------------
    !
    MV = 0.0D0
    !
    DO I = 0, DIM - 1
        !
        MV(:) = MV(:) + M(:, I) * V(I)
        !
    END DO
    !
    END SUBROUTINE MATRIX_VEC_MULT
    !
    !===========================================================================
    !
    SUBROUTINE MATT_X_MAT3(A, B, C, N, M)
    !
    ! Matrix multiplication by transpose. (3x3)
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! A, B: arrays of 3x3 matrices
    ! C: array of 3x3 matrices, c=a^t*b
    ! N: number of grains
    ! M: number of elements
    !
    REAL(RK), INTENT(IN) :: A(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: B(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: C(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    ! I, J, K: Looping indices
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: K
    !
    !----------------------------------------------------------------------
    !
    C = 0.0D0
    !
    DO J = 0, DIMS1
        !
        DO I = 0, DIMS1
            !
            DO K = 0, DIMS1
                !
                C(I, J, :, :) = C(I, J, :, :) + A(K, I, :, :) * B(K, J, :, :)
                !
            END DO
            !
        END DO
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE MATT_X_MAT3
    !
    !===========================================================================
    !
    SUBROUTINE NORM_VEC(NORM, VECTOR, N, M)
    !
    ! Compute 2-norm for array of 5-vectors.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! NORM: Array of vector norms
    ! VECTOR: Array of 5-vectors
    ! N: Number of grains/element (legacy, always 0:0)
    ! M: Number of elements
    !
    REAL(RK), INTENT(OUT) :: NORM(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: VECTOR(0:TVEC1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    ! I: Looping index
    !
    INTEGER :: I
    !
    !----------------------------------------------------------------------
    !
    NORM = 0.0D0
    !
    DO I = 0, TVEC1
        !
        NORM = NORM + VECTOR(I, :, :) * VECTOR(I, :, :)
        !
    END DO
    !
    NORM = DSQRT(NORM)
    !
    RETURN
    !
    END SUBROUTINE NORM_VEC
    !
    !===========================================================================
    !
    SUBROUTINE ROT_MAT_SKEW(C, QR3X3, N, M)
    !
    ! Construct 3x3 matrix acting on skew 3-vectors which satisfy the matrix
    !   operation:
    !
    ! [W]_SAM = [C] * [W]_LAT * [C]'
    !
    ! thus, for vectors:
    !
    ! {W}_SAM = [QR3X3] * {W}_LAT
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! C: Array of orientation matrices
    ! QR3X3: Array of 3x3 matrices which acts on skew vectors
    ! N: Number of grains/element (legacy, always 0:0)
    ! M: Number of elements
    !
    REAL(RK), INTENT(IN)  :: C(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: QR3X3(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    ! Cxx: Intermediary calculations
    !
    REAL(RK) :: C11(0:(N - 1), 0:(M - 1))
    REAL(RK) :: C12(0:(N - 1), 0:(M - 1))
    REAL(RK) :: C13(0:(N - 1), 0:(M - 1))
    REAL(RK) :: C21(0:(N - 1), 0:(M - 1))
    REAL(RK) :: C22(0:(N - 1), 0:(M - 1))
    REAL(RK) :: C23(0:(N - 1), 0:(M - 1))
    REAL(RK) :: C31(0:(N - 1), 0:(M - 1))
    REAL(RK) :: C32(0:(N - 1), 0:(M - 1))
    REAL(RK) :: C33(0:(N - 1), 0:(M - 1))
    !
    !---------------------------------------------------------------------------
    !
    C11 = C(0, 0, :, :)
    C21 = C(1, 0, :, :)
    C31 = C(2, 0, :, :)
    C12 = C(0, 1, :, :)
    C22 = C(1, 1, :, :)
    C32 = C(2, 1, :, :)
    C13 = C(0, 2, :, :)
    C23 = C(1, 2, :, :)
    C33 = C(2, 2, :, :)
    !
    QR3X3(0, 0, :, :) = C22 * C11 - C21 * C12
    QR3X3(0, 1, :, :) = C23 * C11 - C21 * C13
    QR3X3(0, 2, :, :) = C23 * C12 - C22 * C13
    QR3X3(1, 0, :, :) = C32 * C11 - C31 * C12
    QR3X3(1, 1, :, :) = C33 * C11 - C31 * C13
    QR3X3(1, 2, :, :) = C33 * C12 - C32 * C13
    QR3X3(2, 0, :, :) = C32 * C21 - C31 * C22
    QR3X3(2, 1, :, :) = C33 * C21 - C31 * C23
    QR3X3(2, 2, :, :) = C33 * C22 - C32 * C23
    !
    RETURN
    !
    END SUBROUTINE ROT_MAT_SKEW
    !
    !===========================================================================
    !
    SUBROUTINE ROT_MAT_SKEW_SER(C, QR3X3)
    !
    ! Construct 3x3 matrix acting on skew 3-vectors which satisfy the matrix
    !   operation:
    !
    ! [W]_SAM = [C] * [W]_LAT * [C]'
    !
    ! thus, for vectors:
    !
    ! {W}_SAM = [QR3X3] * {W}_LAT
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! C: Orientation matrix
    ! QR3X3: 3x3 matrix which acts on skew vectors
    !
    REAL(RK), INTENT(IN)  :: C(0:DIMS1, 0:DIMS1)
    REAL(RK), INTENT(OUT) :: QR3X3(0:DIMS1, 0:DIMS1)
    !
    ! Locals:
    ! Cxx: Intermediary calculations
    !
    REAL(RK) :: C11
    REAL(RK) :: C12
    REAL(RK) :: C13
    REAL(RK) :: C21
    REAL(RK) :: C22
    REAL(RK) :: C23
    REAL(RK) :: C31
    REAL(RK) :: C32
    REAL(RK) :: C33
    !
    !---------------------------------------------------------------------------
    !
    C11 = C(0, 0)
    C21 = C(1, 0)
    C31 = C(2, 0)
    C12 = C(0, 1)
    C22 = C(1, 1)
    C32 = C(2, 1)
    C13 = C(0, 2)
    C23 = C(1, 2)
    C33 = C(2, 2)
    !
    QR3X3(0, 0) = C22 * C11 - C21 * C12
    QR3X3(0, 1) = C23 * C11 - C21 * C13
    QR3X3(0, 2) = C23 * C12 - C22 * C13
    QR3X3(1, 0) = C32 * C11 - C31 * C12
    QR3X3(1, 1) = C33 * C11 - C31 * C13
    QR3X3(1, 2) = C33 * C12 - C32 * C13
    QR3X3(2, 0) = C32 * C21 - C31 * C22
    QR3X3(2, 1) = C33 * C21 - C31 * C23
    QR3X3(2, 2) = C33 * C22 - C32 * C23
    !
    RETURN
    !
    END SUBROUTINE ROT_MAT_SKEW_SER
    !
    !===========================================================================
    !
    SUBROUTINE ROT_MAT_SYMM(C, QR5X5, N, M)
    !
    ! Convert rotation matrix C to an operator on 5-vectors which satisfy the
    !   matrix operation:
    !
    ! [A]_SAM = [C] * [A]_LAT * [C]'
    !
    ! thus, for vectors:
    !
    ! {A}_SAM = [QR5X5] * {A}_LAT
    !
    ! with: {A}={()/SQR2,SQR32*(),SQR2*(),SQR2*(),SQR2*()}
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! C: Array of rotation matrices
    ! QR5X5: Array of 5x5 rotation matrices
    ! N: Number of grains/element (legacy, always 0:0)
    ! M: Number of elements
    !
    REAL(RK), INTENT(IN)  :: C(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: QR5X5(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    ! SQR3: Local parameter
    ! Cxx: Intermediary calculations
    !
    REAL(RK) :: SQR3
    REAL(RK) :: C11(0:(N - 1), 0:(M - 1))
    REAL(RK) :: C12(0:(N - 1), 0:(M - 1))
    REAL(RK) :: C13(0:(N - 1), 0:(M - 1))
    REAL(RK) :: C21(0:(N - 1), 0:(M - 1))
    REAL(RK) :: C22(0:(N - 1), 0:(M - 1))
    REAL(RK) :: C23(0:(N - 1), 0:(M - 1))
    REAL(RK) :: C31(0:(N - 1), 0:(M - 1))
    REAL(RK) :: C32(0:(N - 1), 0:(M - 1))
    REAL(RK) :: C33(0:(N - 1), 0:(M - 1))
    !
    !----------------------------------------------------------------------
    !
    SQR3 = DSQRT(3.0D0)
    !
    C11 = C(0, 0, :, :)
    C21 = C(1, 0, :, :)
    C31 = C(2, 0, :, :)
    C12 = C(0, 1, :, :)
    C22 = C(1, 1, :, :)
    C32 = C(2, 1, :, :)
    C13 = C(0, 2, :, :)
    C23 = C(1, 2, :, :)
    C33 = C(2, 2, :, :)
    !
    QR5X5(0, 0, :, :)  =  0.5D0 * (C11 * C11 - C12 * C12 - C21 * C21 + C22 * &
        & C22)
    QR5X5(0, 1, :, :)  =  SQR3 / 2.0D0 * (C13 * C13 - C23 * C23)
    QR5X5(0, 2, :, :)  =  C11 * C12 - C21 * C22
    QR5X5(0, 3, :, :)  =  C11 * C13 - C21 * C23
    QR5X5(0, 4, :, :)  =  C12 * C13 - C22 * C23
    QR5X5(1, 0, :, :)  =  SQR3 / 2.0D0 * (C31 * C31 - C32 * C32)
    QR5X5(1, 1, :, :)  =  1.5D0 * C33 * C33 - 0.5D0
    QR5X5(1, 2, :, :)  =  SQR3 * C31 * C32
    QR5X5(1, 3, :, :)  =  SQR3 * C31 * C33
    QR5X5(1, 4, :, :)  =  SQR3 * C32 * C33
    QR5X5(2, 0, :, :)  =  C11 * C21 - C12 * C22
    QR5X5(2, 1, :, :)  =  SQR3 * C13 * C23
    QR5X5(2, 2, :, :)  =  C11 * C22 + C12 * C21
    QR5X5(2, 3, :, :)  =  C11 * C23 + C13 * C21
    QR5X5(2, 4, :, :)  =  C12 * C23 + C13 * C22
    QR5X5(3, 0, :, :)  =  C11 * C31 - C12 * C32
    QR5X5(3, 1, :, :)  =  SQR3 * C13 * C33
    QR5X5(3, 2, :, :)  =  C11 * C32 + C12 * C31
    QR5X5(3, 3, :, :)  =  C11 * C33 + C13 * C31
    QR5X5(3, 4, :, :)  =  C12 * C33 + C13 * C32
    QR5X5(4, 0, :, :)  =  C21 * C31 - C22 * C32
    QR5X5(4, 1, :, :)  =  SQR3 * C23 * C33
    QR5X5(4, 2, :, :)  =  C21 * C32 + C22 * C31
    QR5X5(4, 3, :, :)  =  C21 * C33 + C23 * C31
    QR5X5(4, 4, :, :)  =  C22 * C33 + C23 * C32
    !
    RETURN
    !
    END SUBROUTINE ROT_MAT_SYMM
    !
    !===========================================================================
    !
    SUBROUTINE ROT_MAT_SYMM_SER(C, QR5X5)
    !
    ! Convert rotation matrix C to an operator on 5-vectors which satisfy the
    !   matrix operation:
    !
    ! [A]_SAM = [C] * [A]_LAT * [C]'
    !
    ! thus, for vectors:
    !
    ! {A}_SAM = [QR5X5] * {A}_LAT
    !
    ! with: {A}={()/SQR2,SQR32*(),SQR2*(),SQR2*(),SQR2*()}
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! C: Rotation matrix
    ! QR5X5: 5x5 rotation matrix
    !
    REAL(RK), INTENT(IN) :: C(0:DIMS1, 0:DIMS1)
    REAL(RK), INTENT(OUT) :: QR5X5(0:TVEC1, 0:TVEC1)
    !
    ! Locals:
    ! SQR3: Local parameter
    ! Cxx: Intermediary calculations
    !
    REAL(RK) :: SQR3
    REAL(RK) :: C11
    REAL(RK) :: C12
    REAL(RK) :: C13
    REAL(RK) :: C21
    REAL(RK) :: C22
    REAL(RK) :: C23
    REAL(RK) :: C31
    REAL(RK) :: C32
    REAL(RK) :: C33
    !
    !----------------------------------------------------------------------
    !
    SQR3 = DSQRT(3.0D0)
    !
    C11 = C(0, 0)
    C21 = C(1, 0)
    C31 = C(2, 0)
    C12 = C(0, 1)
    C22 = C(1, 1)
    C32 = C(2, 1)
    C13 = C(0, 2)
    C23 = C(1, 2)
    C33 = C(2, 2)
    !
    QR5X5(0, 0)  =  0.5D0 * (C11 * C11 - C12 * C12 - C21 * C21 + C22 * C22)
    QR5X5(0, 1)  =  SQR3 / 2.0D0 * (C13 * C13 - C23 * C23)
    QR5X5(0, 2)  =  C11 * C12 - C21 * C22
    QR5X5(0, 3)  =  C11 * C13 - C21 * C23
    QR5X5(0, 4)  =  C12 * C13 - C22 * C23
    QR5X5(1, 0)  =  SQR3 / 2.0D0 * (C31 * C31 - C32 * C32)
    QR5X5(1, 1)  =  1.5D0 * C33 * C33 - 0.5D0
    QR5X5(1, 2)  =  SQR3 * C31 * C32
    QR5X5(1, 3)  =  SQR3 * C31 * C33
    QR5X5(1, 4)  =  SQR3 * C32 * C33
    QR5X5(2, 0)  =  C11 * C21 - C12 * C22
    QR5X5(2, 1)  =  SQR3 * C13 * C23
    QR5X5(2, 2)  =  C11 * C22 + C12 * C21
    QR5X5(2, 3)  =  C11 * C23 + C13 * C21
    QR5X5(2, 4)  =  C12 * C23 + C13 * C22
    QR5X5(3, 0)  =  C11 * C31 - C12 * C32
    QR5X5(3, 1)  =  SQR3 * C13 * C33
    QR5X5(3, 2)  =  C11 * C32 + C12 * C31
    QR5X5(3, 3)  =  C11 * C33 + C13 * C31
    QR5X5(3, 4)  =  C12 * C33 + C13 * C32
    QR5X5(4, 0)  =  C21 * C31 - C22 * C32
    QR5X5(4, 1)  =  SQR3 * C23 * C33
    QR5X5(4, 2)  =  C21 * C32 + C22 * C31
    QR5X5(4, 3)  =  C21 * C33 + C23 * C31
    QR5X5(4, 4)  =  C22 * C33 + C23 * C32
    !
    RETURN
    !
    END SUBROUTINE ROT_MAT_SYMM_SER
    !
    !===========================================================================
    !
    SUBROUTINE SOLVE_LIN_SYS_3(MAT, VEC, SOL)
    !
    ! Solve a linear system of three equations. Used for triaxial loading.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    REAL(RK), INTENT(IN) :: MAT(3,3)
    REAL(RK), INTENT(IN) :: VEC(3)
    REAL(RK), INTENT(OUT) :: SOL(3)
    !
    ! Locals:
    !
    REAL(RK) :: INV(3,3)
    REAL(RK) :: A, B, C, D, E, F, G, H, K
    REAL(RK) :: DET, MATNORM, INVNORM, COND
    !
    !---------------------------------------------------------------------------
    !
    A = MAT(1, 1)
    B = MAT(1, 2)
    C = MAT(1, 3)
    D = MAT(2, 1)
    E = MAT(2, 2)
    F = MAT(2, 3)
    G = MAT(3, 1)
    H = MAT(3, 2)
    K = MAT(3, 3)
    !
    DET = A * (E * K - F * H) - B * (D * K - F * G) + C * (D * H - E * G)
    !
    INV(1, 1) = (E * K - F * H) / DET
    INV(1, 2) = (C * H - B * K) / DET
    INV(1, 3) = (B * F - C * E) / DET
    INV(2, 1) = (F * G - D * K) / DET
    INV(2, 2) = (A * K - C * G) / DET
    INV(2, 3) = (C * D - A * F) / DET
    INV(3, 1) = (D * H - E * G) / DET
    INV(3, 2) = (B * G - A * H) / DET
    INV(3, 3) = (A * E - B * D) / DET
    !
    SOL = MATMUL(INV, VEC)
    !
    ! Find conditioning number
    !
    MATNORM = MAX(MAT(1, 1) + MAT(1, 2) + MAT(1, 3), &
         & MAT(2, 1) + MAT(2, 2) + MAT(2, 3), &
         & MAT(3, 1) + MAT(3, 2) + MAT(3, 3))
    INVNORM = MAX(INV(1, 1) + INV(1, 2) + INV(1, 3), &
         & INV(2, 1) + INV(2, 2) + INV(2, 3), &
         & INV(3, 1) + INV(3, 2) + INV(3, 3))
    COND = MATNORM * INVNORM
    !
    IF (COND .GT. 1.0D3) THEN
        !
        CALL PAR_QUIT('Error  :     > Matrix is poorly conditioned.')
        !
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE SOLVE_LIN_SYS_3
    !
    !===========================================================================
    !
    SUBROUTINE SPARSE_MATVEC_EBE(RES, SOL, TEMP1, TEMP2, GSTIF, BCS, NNPE, &
        & NSUB, NSUP, ESUB, ESUP, DTRACE, NP)
    !
    !  Matrix vector multiply, element by element.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! RES:
    ! SOL:
    ! TEMP1:
    ! TEMP2:
    ! GSTIF:
    ! BCS:
    ! NNPE:
    ! NSUB:
    ! NSUP:
    ! ESUB:
    ! ESUP:
    ! DTRACE:
    ! NP:
    !
    REAL(RK) :: RES(NSUB:NSUP)
    REAL(RK) :: SOL(NSUB:NSUP)
    REAL(RK) :: TEMP1(0:(NNPE - 1), ESUB:ESUP)
    REAL(RK) :: TEMP2(0:(NNPE - 1), ESUB:ESUP)
    REAL(RK) :: GSTIF(0:(NNPE - 1), 0:(NNPE - 1), ESUB:ESUP)
    LOGICAL :: BCS(NSUB:NSUP)
    INTEGER :: NNPE
    INTEGER :: NSUB
    INTEGER :: NSUP
    INTEGER :: ESUB
    INTEGER :: ESUP
    TYPE(TRACE) :: DTRACE
    INTEGER :: NP(0:(NNPE-1),ESUB:ESUP)
    !
    !  Locals:
    !
    INTEGER :: IER
    INTEGER :: I
    INTEGER :: J
    INTEGER :: IDUMMY
    !
    !---------------------------------------------------------------------------
    !
    CALL PART_GATHER(TEMP1, SOL, NP, DTRACE)
    !
    TEMP2 = 0.0D0
    !
    CALL GEN_MATRIX_VECTOR_MULT(TEMP2, GSTIF, TEMP1, IDUMMY, IDUMMY, IDUMMY, &
        & IDUMMY, IER)
    !
    RES = 0.0D0
    !
    CALL PART_SCATTER(RES, TEMP2, NP, .FALSE., DTRACE)
    !
    WHERE (BCS)
        !
        RES = 0.0D0
        !
    END WHERE
    !
    RETURN
    !
    END SUBROUTINE SPARSE_MATVEC_EBE
    !
    !===========================================================================
    !
    SUBROUTINE SYMM_VGR(D, DKK, VGRAD, M)
    !
    ! Compute symmetric part of an array of velocity gradients.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! D: Deviatoric part of velocity gradients
    ! DKK: Trace of velocity gradients
    ! VGRAD: Array of velocity gradients
    ! M: Number of elements
    !
    REAL(RK), INTENT(OUT) :: D(0:DIMS1, 0:DIMS1, 0:(M - 1))
    REAL(RK), INTENT(OUT) :: DKK(0:(M - 1))
    REAL(RK), INTENT(IN) :: VGRAD(0:DIMS1, 0:DIMS1, 0:(m -1))
    INTEGER, INTENT(IN) :: M
    !
    !---------------------------------------------------------------------------
    !
    D(0, 0, :) = VGRAD(0, 0, :)
    D(1, 1, :) = VGRAD(1, 1, :)
    D(2, 2, :) = VGRAD(2, 2, :)
    D(1, 0, :) = 0.5D0 * (VGRAD(1, 0, :) + VGRAD(0, 1, :))
    D(2, 0, :) = 0.5D0 * (VGRAD(2, 0, :) + VGRAD(0, 2, :))
    D(2, 1, :) = 0.5D0 * (VGRAD(2, 1, :) + VGRAD(1, 2, :))
    D(0, 1, :) = D(1, 0, :)
    D(0, 2, :) = D(2, 0, :)
    D(1, 2, :) = D(2, 1, :)
    !
    DKK = D(0, 0, :) + D(1, 1, :) + D(2, 2, :)
    !
    D(0, 0, :) = D(0, 0, :) - DKK / 3.0D0
    D(1, 1, :) = D(1, 1, :) - DKK / 3.0D0
    D(2, 2, :) = D(2, 2, :) - DKK / 3.0D0
    !
    RETURN
    !
    END SUBROUTINE SYMM_VGR
    !
    !===========================================================================
    !
    SUBROUTINE SYMM_VGR_SER(D, DKK, VGRAD)
    !
    ! Compute symmetric part of velocity gradients.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! D: deviatoric part of velocity gradients
    ! DKK: trace of velocity gradients
    ! VGRAD: array of velocity gradients
    ! M: number of elements
    !
    REAL(RK), INTENT(OUT) :: D(0:DIMS1, 0:DIMS1)
    REAL(RK), INTENT(OUT) :: DKK
    REAL(RK), INTENT(IN) :: VGRAD(0:DIMS1, 0:DIMS1)
    !
    !---------------------------------------------------------------------------
    !
    D(0, 0) = VGRAD(0, 0)
    D(1, 1) = VGRAD(1, 1)
    D(2, 2) = VGRAD(2, 2)
    D(1, 0) = 0.5D0 * (VGRAD(1, 0) + VGRAD(0, 1))
    D(2, 0) = 0.5D0 * (VGRAD(2, 0) + VGRAD(0, 2))
    D(2, 1) = 0.5D0 * (VGRAD(2, 1) + VGRAD(1, 2))
    D(0, 1) = D(1, 0)
    D(0, 2) = D(2, 0)
    D(1, 2) = D(2, 1)
    !
    DKK = D(0, 0) + D(1, 1) + D(2, 2)
    !
    D(0, 0) = D(0, 0) - DKK / 3.0D0
    D(1, 1) = D(1, 1) - DKK / 3.0D0
    D(2, 2) = D(2, 2) - DKK / 3.0D0
    !
    RETURN
    !
    END SUBROUTINE SYMM_VGR_SER
    !
    !===========================================================================
    !
    SUBROUTINE SKEW_VGR(W, VGRAD, M)
    !
    ! Compute skew part of an array of velocity gradients.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! W: Array of skew matrices (output)
    ! VGRAD: Array of velocity gradients
    ! M: Number of elements
    !
    REAL(RK), INTENT(IN) :: VGRAD(0:DIMS1, 0:DIMS1, 0:(M - 1))
    REAL(RK), INTENT(OUT) :: W(0:DIMS1, 0:DIMS1, 0:(M - 1))
    INTEGER, INTENT(IN) :: M
    !
    !---------------------------------------------------------------------------
    !
    W(:, :, :) = 0.0D0
    !
    W(0, 1, :) = 0.5D0 * (VGRAD(0, 1, :) - VGRAD(1, 0, :))
    W(0, 2, :) = 0.5D0 * (VGRAD(0, 2, :) - VGRAD(2, 0, :))
    W(1, 2, :) = 0.5D0 * (VGRAD(1, 2, :) - VGRAD(2, 1, :))
    W(1, 0, :) = -W(0, 1, :)
    W(2, 0, :) = -W(0, 2, :)
    W(2, 1, :) = -W(1, 2, :)
    !
    RETURN
    !
    END SUBROUTINE SKEW_VGR
    !
    !===========================================================================
    !
    SUBROUTINE TENSOR3DCOMPOSE(MAT, DEV, SKW, SPH)
    !
    ! Form  matrix from deviatoric, skew, or spherical parts. If no parts are
    !   passed, the matix is zeroed.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! MAT: Resulting array of 3x3 matrices
    ! DEV: Arry of 5-vec representing symmetric, traceless portion
    ! SKW: Array of 3-vec representing axial vectors of the skew pportion
    ! SPH: Array of scalars representing one third of the trace
    !
    REAL(RK), INTENT(OUT) :: MAT(:, :, :)
    REAL(RK), INTENT(IN), OPTIONAL :: DEV(:, :)
    REAL(RK), INTENT(IN), OPTIONAL :: SKW(:, :)
    REAL(RK), INTENT(IN), OPTIONAL :: SPH(:)
    !
    !---------------------------------------------------------------------------
    !
    MAT = 0.0D0
    !
    IF (PRESENT(DEV)) THEN
        !
        MAT(1, 1, :) = -DEV(1, :) * SQ2_I - DEV(2, :) * SQ6_I
        MAT(2, 2, :) = DEV(1, :) * SQ2_I - DEV(2, :) * SQ6_I
        MAT(3, 3, :) = DEV(2, :) * TWOSQ6_I
        !
        MAT(2, 3, :) = DEV(3, :) * SQ2_I
        MAT(3, 2, :) = MAT(2, 3, :)
        !
        MAT(1, 3, :) = DEV(4, :) * SQ2_I
        MAT(3, 1, :) = MAT(1, 3, :)
        !
        MAT(1, 2, :) = DEV(5, :) * SQ2_I
        MAT(2, 1, :) = MAT(1, 2, :)
        !
    END IF
    !
    IF (PRESENT(SKW)) THEN
        !
        MAT(2, 3, :) = MAT(2, 3, :) - SKW(1, :)
        MAT(3, 2, :) = MAT(3, 2, :) + SKW(1, :)
        !
        MAT(1, 3, :) = MAT(1, 3, :) + SKW(2, :)
        MAT(3, 1, :) = MAT(3, 1, :) - SKW(2, :)
        !
        MAT(1, 2, :) = MAT(1, 2, :) - SKW(3, :)
        MAT(2, 1, :) = MAT(2, 1, :) + SKW(3, :)
        !
    END IF
    !
    IF (PRESENT(SPH)) THEN
        !
        MAT(1, 1, :) = MAT(1, 1, :) + SPH
        MAT(2, 2, :) = MAT(2, 2, :) + SPH
        MAT(3, 3, :) = MAT(3, 3, :) + SPH
        !
    END IF
    !
    END SUBROUTINE TENSOR3DCOMPOSE
    !
    !===========================================================================
    !
    SUBROUTINE TENSOR3DDECOMPOSE(MAT, DEV, SKW, SPH, DECOMP)
    !
    ! Decompose matrix into deviatoric, skew, or spherical parts.
    !
    ! Note that the three components of the decomposition are orthogonal, but
    !   the underlying basis is not orthonormal due to scaling.
    !
    ! Note: dev(2) = sqrt(3/2)*mat(3,3), when mat is deviatoric.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! MAT: Array of 3x3 matrices
    ! DEV: Array of 5-vec representing symmetric portion (MPSIM convention)
    ! SKW: Array of 3-vec representing skew portion (MPSIM convention)
    ! SPH: Array of scalars representing one third of the trace
    ! DECOMP: Flag indicating decomposition convention
    !
    REAL(RK), INTENT(IN) :: MAT(:, :, :)
    REAL(RK), INTENT(OUT), OPTIONAL :: DEV(:, :)
    REAL(RK), INTENT(OUT), OPTIONAL :: SKW(:, :)
    REAL(RK), INTENT(OUT), OPTIONAL :: SPH(:)
    INTEGER, INTENT(IN), OPTIONAL :: DECOMP
    !
    ! Locals:
    !
    INTEGER ::  DECOMP_CONV
    !
    !---------------------------------------------------------------------------
    !
    DECOMP_CONV = DECOMP_DFLT
    !
    IF (PRESENT(DECOMP)) THEN
        !
        DECOMP_CONV = DECOMP
        !
    END IF
    !
    IF (PRESENT(DEV)) THEN
        !
        SELECT CASE(DECOMP_CONV)
        !
        CASE (DECOMP_MPSIM)
            !
            DEV(1, :) = (MAT(2, 2, :) - MAT(1, 1, :)) * SQ2_I
            DEV(2, :) = (MAT(3, 3, :) + MAT(3, 3, :) - MAT(1, 1, :) - &
                & MAT(2, 2, :)) * SQ6_I
            !
            DEV(3, :) = SQ2_I * (MAT(2, 3, :) + MAT(3, 2, :))
            DEV(4, :) = SQ2_I * (MAT(1, 3, :) + MAT(3, 1, :))
            DEV(5, :) = SQ2_I * (MAT(1, 2, :) + MAT(2, 1, :))
            !
        CASE (DECOMP_FEMEVPS)
            !
            DEV(1,:) = ( MAT(1,1,:) - MAT(2,2,:) )*SQ2_I
            DEV(2,:) = ( MAT(3,3,:) + MAT(3,3,:) - MAT(1,1,:) - MAT(2,2,:) ) * &
                & SQ6_I
            !
            DEV(3, :) = SQ2_I * (MAT(1, 2, :) + MAT(2, 1, :))
            DEV(4, :) = SQ2_I * (MAT(1, 3, :) + MAT(3, 1, :))
            DEV(5, :) = SQ2_I * (MAT(2, 3, :) + MAT(3, 2, :))
            !
        CASE DEFAULT
        !
        ! Return error status
        !
        END SELECT
        !
    END IF
    !
    IF (PRESENT(SKW)) THEN
        !
        SELECT CASE(DECOMP_CONV)
        !
        CASE (DECOMP_MPSIM)
            !
            SKW(1, :) = 0.5D0 * (MAT(3, 2, :) - MAT(2, 3, :))
            SKW(2, :) = 0.5D0 * (MAT(1, 3, :) - MAT(3, 1, :))
            SKW(3, :) = 0.5D0 * (MAT(2, 1, :) - MAT(1, 2, :))
            !
        CASE (DECOMP_FEMEVPS)
            !
            SKW(1, :) = 0.5D0 * (MAT(2, 1, :) - MAT(1, 2, :))
            SKW(2, :) = -0.5D0 * (MAT(1, 3, :) - MAT(3, 1, :))
            SKW(3, :) = 0.5D0 * (MAT(3, 2, :) - MAT(2, 3, :))
            !
        CASE DEFAULT
            !
            ! Return error status
            !
        END SELECT
        !
    END IF
    !
    IF (PRESENT(SPH)) THEN
        !
        SPH = (MAT(1, 1, :) + MAT(2, 2, :) + MAT(3, 3, :)) * (1.0D0 / 3.0D0)
        !
    END IF
    !
    END SUBROUTINE TENSOR3DDECOMPOSE
    !
    !===========================================================================
    !
    SUBROUTINE VEC5_VEC6(VEC, MAT)
    !
    ! Convert 5-vector to 6-vector of symmetric matrix.
    !
    ! Note: Returns the upper triangle in the order of: 11 22 33 23 13 12
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! VEC: 5-vector of dev. part of 2nd-order tensor (inner product preserved)
    ! MAT: 6-vector of upper triangle of symmetric tensor (no IP preservation)
    !
    REAL(RK), INTENT(IN) :: VEC(0:TVEC1)
    REAL(RK), INTENT(OUT) :: MAT(0:TVEC)
    !
    ! Locals:
    ! SQR2, SQR23: Local parameters
    !
    REAL(RK) :: SQR2
    REAL(RK) :: SQR23
    !
    !---------------------------------------------------------------------------
    !
    SQR2 = DSQRT(2.0D0)
    SQR23 = DSQRT(2.0D0 / 3.0D0)
    !
    ! Diagonal: 11 22 33
    MAT(0) = 0.5D0 * (SQR2 * VEC(0) - SQR23 * VEC(1))
    MAT(1) = -0.5D0 * (SQR2 * VEC(0) + SQR23 * VEC(1))
    MAT(2) = VEC(1) * SQR23
    !
    ! Off-diagonal: 23 13 12
    MAT(3) = VEC(4) / SQR2
    MAT(4) = VEC(3) / SQR2
    MAT(5) = VEC(2) / SQR2
    !
    RETURN
    !
    END SUBROUTINE VEC5_VEC6
    !
    !===========================================================================
    !
    SUBROUTINE VEC_D_VEC5(D, B, C, N, M)
    !
    ! Computes c = d*b where c and b are 5-vectors and d is a diagonal 5x5
    !   matrix represented in compact form
    !
    !----------------------------------------------------------------------
    !
    ! Arguments:
    ! D: A fixed diagonal matrix (in 5-vector form)
    ! B: Input matrix (in 5-vector form)
    ! C: Output matrix (d*b) (in 5-vector form)
    ! N: Number of grains/element (legacy, always 0:0)
    ! M: Number of elements
    !
    REAL(RK), INTENT(IN) :: D(0:TVEC1)
    REAL(RK), INTENT(IN) :: B(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: C(0:TVEC1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    ! K: Looping index
    !
    INTEGER :: K
    !
    !----------------------------------------------------------------------
    !
    DO K = 0, TVEC1
        !
        C(K, :, :) = D(K) * B(K, :, :)
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE VEC_D_VEC5
    !
    !===========================================================================
    !
    SUBROUTINE VEC_MAT_SKEW(VEC, MAT, M)
    !
    ! Convert array of 3-vectors to array of skew matrices.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! VEC: Array of 3-vectors
    ! MAT: Array of skew matrices
    ! M: Number of elements
    !
    REAL(RK), INTENT(IN) :: VEC(0:DIMS1, 0:(M - 1))
    REAL(RK), INTENT(OUT) :: MAT(0:DIMS1, 0:DIMS1, 0:(M - 1))
    INTEGER, INTENT(IN) :: M
    !
    !---------------------------------------------------------------------------
    !
    MAT(0, 0, :) = 0.0D0
    MAT(1, 1, :) = 0.0D0
    MAT(2, 2, :) = 0.0D0
    !
    MAT(1, 0, :) = VEC(0, :)
    MAT(2, 0, :) = VEC(1, :)
    MAT(2, 1, :) = VEC(2, :)

    MAT(0, 1, :) = -VEC(0, :)
    MAT(0, 2, :) = -VEC(1, :)
    MAT(1, 2, :) = -VEC(2, :)
    !
    RETURN
    !
    END SUBROUTINE VEC_MAT_SKEW
    !
    !===========================================================================
    !
    SUBROUTINE VEC_MAT_SKEW_GRN(VEC, MAT, N, M)
    !
    ! Convert 3-vectors to skew matrices across grains and elements.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! VEC: Array of 3-vectors
    ! MAT: Array of skew matrices
    ! N: Number of grains/element (legacy, always 0:0)
    ! M: Number of elements
    !
    
    REAL(RK), INTENT(IN) :: VEC(0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: MAT(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    !----------------------------------------------------------------------
    !
    MAT(0, 0, :, :) = 0.0D0
    MAT(1, 1, :, :) = 0.0D0
    MAT(2, 2, :, :) = 0.0D0
    !
    MAT(1, 0, :, :) = VEC(0, :, :)
    MAT(2, 0, :, :) = VEC(1, :, :)
    MAT(2, 1, :, :) = VEC(2, :, :)
    !
    MAT(0, 1, :, :) = -VEC(0, :, :)
    MAT(0, 2, :, :) = -VEC(1, :, :)
    MAT(1, 2, :, :) = -VEC(2, :, :)
    !
    RETURN
    !
    END SUBROUTINE VEC_MAT_SKEW_GRN
    !
    !===========================================================================
    !
    SUBROUTINE VEC_MAT_SYMM(VEC, MAT, M)
    !
    ! Convert 5-vector to symmetric matrix: {5} --> [3x3]sym
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! VEC: Array of 5-vectors
    ! MAT: Array of 3x3 symmetric matrices
    ! M: number of elements
    !
    REAL(RK), INTENT(IN) :: VEC(0:TVEC1, 0:(M - 1))
    REAL(RK), INTENT(OUT) :: MAT(0:DIMS1, 0:DIMS1, 0:(M - 1))
    INTEGER , INTENT(IN) :: M
    !
    ! Locals:
    ! SQR2, SQR23: Local parameters
    !
    REAL(RK) :: SQR2
    REAL(RK) :: SQR23
    !
    !----------------------------------------------------------------------
    !
    SQR2 = DSQRT(2.0D0)
    SQR23 = DSQRT(2.0D0 / 3.0D0)
    !
    MAT(0, 0, :) = 0.5D0 * (SQR2 * VEC(0, :) - SQR23 * VEC(1, :))
    MAT(1, 1, :) = -0.5D0 * (SQR2 * VEC(0, :) + SQR23 * VEC(1, :))
    MAT(2, 2, :) = VEC(1, :) * SQR23
    MAT(1, 0, :) = VEC(2, :) / SQR2
    MAT(2, 0, :) = VEC(3, :) / SQR2
    MAT(2, 1, :) = VEC(4, :) / SQR2
    !
    MAT(0, 1, :) = MAT(1, 0, :)
    MAT(0, 2, :) = MAT(2, 0, :)
    MAT(1, 2, :) = MAT(2, 1, :)
    !
    RETURN
    !
    END SUBROUTINE VEC_MAT_SYMM
    !
    !===========================================================================
    !
    SUBROUTINE VEC_MAT_SYMM_GRN(VEC, MAT, N, M)
    !
    ! Convert 5-vector to symmetric matrix: {5} --> [3x3]sym
    !
    !----------------------------------------------------------------------
    !
    ! Arguments:
    ! VEC: Array of 5-vectors
    ! MAT: Array of 3x3 symmetric matrices
    ! N: Number of grains/element (legacy, always 0:0)
    ! M: Number of elements
    !
    REAL(RK), INTENT(IN)  :: VEC(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: MAT(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    INTEGER , INTENT(IN)  :: N
    INTEGER , INTENT(IN)  :: M
    !
    ! Locals:
    ! SQR2, SQR23: Local parameters
    !
    REAL(RK) :: SQR2
    REAL(RK) :: SQR23
    !
    !----------------------------------------------------------------------
    !
    SQR2 = DSQRT(2.0D0)
    SQR23 = DSQRT(2.0D0 / 3.0D0)
    !
    MAT(0, 0, :, :) = 0.5D0 * (SQR2 * VEC(0, :, :) - SQR23 * VEC(1, :, :))
    MAT(1, 1, :, :) = -0.5D0 * (SQR2 * VEC(0, :, :) + SQR23 * VEC(1, :, :))
    MAT(2, 2, :, :) = VEC(1, :, :) * SQR23
    MAT(1, 0, :, :) = VEC(2, :, :) / SQR2
    MAT(2, 0, :, :) = VEC(3, :, :) / SQR2
    MAT(2, 1, :, :) = VEC(4, :, :) / SQR2
    !
    MAT(0, 1, :, :) = MAT(1, 0, :, :)
    MAT(0, 2, :, :) = MAT(2, 0, :, :)
    MAT(1, 2, :, :) = MAT(2, 1, :, :)
    !
    RETURN
    !
    END SUBROUTINE VEC_MAT_SYMM_GRN
    !
    !===========================================================================
    !
    SUBROUTINE VEC_MAT_SYMM_SER(VEC, MAT)
    !
    ! Convert 5-vector to symmetric matrix: {5} --> [3x3]sym
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! VEC: 5-vector
    ! MAT: 3x3 symmetric matrix
    !
    REAL(RK), INTENT(IN) :: VEC(0:TVEC1)
    REAL(RK), INTENT(OUT) :: MAT(0:DIMS1, 0:DIMS1)
    !
    ! Locals:
    ! SQR2, SQR23: Local parameters
    !
    REAL(RK) :: SQR2
    REAL(RK) :: SQR23
    !
    !---------------------------------------------------------------------------
    !
    SQR2 = DSQRT(2.0D0)
    SQR23 = DSQRT(2.0D0 / 3.0D0)
    !
    MAT(0, 0) = 0.5D0 * (SQR2 * VEC(0) - SQR23 * VEC(1))
    MAT(1, 1) = -0.5D0 * (SQR2 * VEC(0) + SQR23 * VEC(1))
    MAT(2, 2) = VEC(1) * SQR23
    MAT(1, 0) = VEC(2) / SQR2
    MAT(2, 0) = VEC(3) / SQR2
    MAT(2, 1) = VEC(4) / SQR2
    !
    MAT(0, 1) = MAT(1, 0)
    MAT(0, 2) = MAT(2, 0)
    MAT(1, 2) = MAT(2, 1)
    !
    RETURN
    !
    END SUBROUTINE VEC_MAT_SYMM_SER
    !
    !===========================================================================
    !
END MODULE MATRIX_OPERATIONS_MOD
