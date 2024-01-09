! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module matrix_operations_mod

! This module contains a variety of matrix operations and other utility
!   functions.

! Contains subroutines:
! find_indices: Find the members of vector equal to target, return val and ind.
! gen_matrix_vector_mult: Matrix times vector
! gen_matrix_mult: Matrix times matrix
! invert5x5: Inverts an array of 5x5 matrices.
! lattice_deform: Change coordinates from sample to lattice (for 5-vectors).
! lattice_spin: Convert a skew 3-vector from sample to lattice coordinates.
! mat_vec_skew: Convert array of skew matrices to array of 3-vectors.
! mat_vec_symm: Convert array of 3x3 symmetric matrices to array of 5-vectors.
! mat_x_mat5: Matrix multiplication for arrays of 5x5 matrices.
! mat_x_matt3: Matrix multiplication by transpose. (3x3)
! mat_x_matt5: Matrix multiplication by transpose for arrays of 5x5 matrices.
! mat_x_mats3: Multiply array of 3x3 matrices by a fixed 3x3 matrix, tranposed.
! mat_x_vec5: Multiply array of matrices times array of vectors. (5 dim)
! matrix_vec_mult: Matrix-vector multiplication routine (rc, 2017)
! vec3_norm: Compute 2-norm for array of 3-vectors.
! vec5_norm: Compute 2-norm for array of 5-vectors.
! rotmat_skew: Construct 3x3 rotation matrix acting on skew-matrix 3-vectors.
! solve_lin_sys_3: Solve a linear system of three equations. Used for triaxial
!   loading.
! sparse_matvec_ebe: Matrix times vector, element by element
! velgrad_sympart: Compute symmetric part of an array of vel gradients.
! velgrad_skewpart: Compute skew part of an array of vel gradients.
! strain_equiv_3x3: Find the equivalent (scalar) value of strain tensors, 3x3.
! stress_equiv_3x3: Find the equivalent (scalar) value for stress tensor, 3x3.
! tensor3dcompose: Form matrix from deviatoric, skew, or spherical parts. If no
!   parts are passed, the matix is zeroed.
! tensor3ddecompose: Decompose matrix into deviatoric and spherical parts.
! vec5_vec6: Convert 5-vector to 6-vector of symmetric matrix.
! vec6_crystosam: Convert 6-vector from crystal to 6-vector sample basis.
! vec_d_vec5: Multiply diagonal matrix times array of vectors. (5 dim)
! vec_mat_skew: Convert array of 3-vectors to array of skew matrices.
! vec_mat_symm: Convert 5-vector to symmetric matrix.

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use quadrature_mod
  use shape_3d_mod
  use gather_scatter_mod

  implicit none

! Constants (all private)

  real(rk), parameter, private :: sq2_i = 1.0d0/dsqrt(2.0d0)
  real(rk), parameter, private :: sq6_i = 1.0d0/dsqrt(6.0d0)
  real(rk), parameter, private :: twosq6_i = 2.0d0*dsqrt(6.0d0)

  public

contains

  !===========================================================================

  subroutine determinant(tensor, determ)

    ! Calculates determinants of an array of 3x3 matrices.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! tensor: Array of 3x3 matrices
    ! determ: Array of scalars (determinant of tensor)
    ! m: Number of elements

    real(rk), intent(in) :: tensor(:, :, :)
    real(rk), intent(out) :: determ(:)

    ! Locals:
    ! DETxx: Intermediary calculations

    real(rk), allocatable :: det11(:)
    real(rk), allocatable :: det12(:)
    real(rk), allocatable :: det13(:)
    real(rk), allocatable :: det21(:)
    real(rk), allocatable :: det22(:)
    real(rk), allocatable :: det23(:)

    integer :: m

    m = size(determ)

    allocate(det11(m))
    allocate(det12(m))
    allocate(det13(m))
    allocate(det21(m))
    allocate(det22(m))
    allocate(det23(m))
    !---------------------------------------------------------------------------

    det11 = tensor(1, 1, :)*tensor(2, 2, :)*tensor(3, 3, :)
    det12 = tensor(1, 2, :)*tensor(2, 3, :)*tensor(3, 1, :)
    det13 = tensor(1, 3, :)*tensor(2, 1, :)*tensor(3, 2, :)
    det21 = tensor(1, 1, :)*tensor(2, 3, :)*tensor(3, 2, :)
    det22 = tensor(1, 2, :)*tensor(2, 1, :)*tensor(3, 3, :)
    det23 = tensor(1, 3, :)*tensor(2, 2, :)*tensor(3, 1, :)

    determ = det11 + det12 + det13 - det21 - det22 - det23

    return

  end subroutine determinant

  !===========================================================================

  ! subroutine find_indices(num_ind, target, vector, indices, offset)
  subroutine find_indices(vector, target, indices, num_ind, offset)

    ! Find the members of vector (integer) which are equal to target (integer),
    !   and return indices and num_ind

    !---------------------------------------------------------------------------

    ! Arguments:

    integer, intent(in) :: vector(:)
    integer, intent(in) :: target
    integer, pointer, intent(out) :: indices(:)
    integer, intent(out) :: num_ind
    integer, optional, intent(in) :: offset

    ! Locals:

    integer :: k
    integer :: lb
    integer :: ub
    integer :: i

    !---------------------------------------------------------------------------

    ! Tito Marin (legacy):

    nullify (indices)

    ! See if indices has been allocated. If so, deallocate.

    if (associated(indices)) then
      deallocate (indices)
    end if

    ! Find how many indices there are.

    num_ind = count(vector .eq. target)

    ! Allocate indices to be num_ind long.

    allocate (indices(num_ind))

    ! Find which indices of vector equal target, and save them in indices.

    k = 0
    lb = lbound(vector, 1) ! Always equals 1
    ub = ubound(vector, 1)

    do i = lb, ub
      if (vector(i) .eq. target) then
        k = k + 1
        indices(k) = i
      end if
    end do

    if (present(offset)) then
      indices = indices + offset
    end if

    return

  end subroutine find_indices

  !===========================================================================

  subroutine gen_matrix_vector_mult(a, x, y)

    !  Array of matrices times array of vectors: y(i) = a(i)*x(i)

    ! Note that the arrays x,y,a are all assumed shape arrays and that they are
    !   therefore dimensioned 1:p, 1:q, etc., even though they may be
    !   dimensioned differently in the calling routine.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! y:
    ! a:
    ! x:

    real(rk), intent(in) :: a(:, :, :)
    real(rk), intent(in) :: x(:, :)
    real(rk), intent(out) :: y(:, :)

    ! Locals:

    integer :: i
    integer :: j
    integer :: ad2
    integer :: ad3

    !---------------------------------------------------------------------------

    ad2 = ubound(a, 2)
    ad3 = ubound(a, 3)

    y = 0.0d0
    do i = 1, ad3
      do j = 1, ad2
        y(:, i) = y(:, i) + a(:, j, i)*x(j, i)
      end do
    end do

    return

  end subroutine gen_matrix_vector_mult

  !===========================================================================

  subroutine gen_matrix_mult(a, b, c)

    ! Accumulate Array of matrices times array of matrices: c(i) = a(i) * b(i)

    ! Note: c is not zeroed in this subroutine, so the values are accumulated;
    !   in all the calling routines, however, the c array is zeroed before
    !   calling this routine.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! c: Array of result matrices (c=a*b)
    ! a, b: Array of input matrices

    real(rk) :: c(:, :, :)
    real(rk) :: a(:, :, :)
    real(rk) :: b(:, :, :)

    ! Locals:

    integer :: i
    integer :: j
    integer :: k
    integer :: lda
    integer :: ldb
    integer :: ldc
    integer :: lta
    integer :: ltb
    integer :: n

    !---------------------------------------------------------------------------

    c = 0.0d0

    lda = ubound(a, 1)
    ldb = ubound(b, 1)
    ldc = ubound(c, 1)
    lta = ubound(a, 2)
    ltb = ubound(b, 2)
    n = ubound(a, 3)

    do i = 1, n
      do j = 1, ltb
        do k = 1, lta
!               do i = 1,ldc
          c(:, j, i) = c(:, j, i) + a(:, k, i)*b(k, j, i)
!               end do
        end do
      end do
    end do

  end subroutine gen_matrix_mult

  !===========================================================================

  subroutine invert5x5(a)

    ! Invert array of 5x5 matrices.

    !---------------------------------------------------------------------------

    ! Arguments:

    ! a: Array of matrices on input and array of inverses on output
    ! m: Number of elements

    real(rk), intent(inout) :: a(:, :, :)

    ! Locals:
    ! i, j: Looping indices
    ! Axx, AIxx: Intermediary calculations

    integer :: i
    integer :: j
    real(rk), allocatable :: a11(:)
    real(rk), allocatable :: a21(:)
    real(rk), allocatable :: a22(:)
    real(rk), allocatable :: a31(:)
    real(rk), allocatable :: a32(:)
    real(rk), allocatable :: a33(:)
    real(rk), allocatable :: a41(:)
    real(rk), allocatable :: a42(:)
    real(rk), allocatable :: a43(:)
    real(rk), allocatable :: a44(:)
    real(rk), allocatable :: a51(:)
    real(rk), allocatable :: a52(:)
    real(rk), allocatable :: a53(:)
    real(rk), allocatable :: a54(:)
    real(rk), allocatable :: a55(:)
    real(rk), allocatable :: ai11(:)
    real(rk), allocatable :: ai21(:)
    real(rk), allocatable :: ai22(:)
    real(rk), allocatable :: ai31(:)
    real(rk), allocatable :: ai32(:)
    real(rk), allocatable :: ai33(:)
    real(rk), allocatable :: ai41(:)
    real(rk), allocatable :: ai42(:)
    real(rk), allocatable :: ai43(:)
    real(rk), allocatable :: ai44(:)
    real(rk), allocatable :: ai51(:)
    real(rk), allocatable :: ai52(:)
    real(rk), allocatable :: ai53(:)
    real(rk), allocatable :: ai54(:)
    real(rk), allocatable :: ai55(:)

    integer :: m

    m = size(a, 3)

    allocate (a11(m))
    allocate (a21(m))
    allocate (a22(m))
    allocate (a31(m))
    allocate (a32(m))
    allocate (a33(m))
    allocate (a41(m))
    allocate (a42(m))
    allocate (a43(m))
    allocate (a44(m))
    allocate (a51(m))
    allocate (a52(m))
    allocate (a53(m))
    allocate (a54(m))
    allocate (a55(m))
    allocate (ai11(m))
    allocate (ai21(m))
    allocate (ai22(m))
    allocate (ai31(m))
    allocate (ai32(m))
    allocate (ai33(m))
    allocate (ai41(m))
    allocate (ai42(m))
    allocate (ai43(m))
    allocate (ai44(m))
    allocate (ai51(m))
    allocate (ai52(m))
    allocate (ai53(m))
    allocate (ai54(m))
    allocate (ai55(m))

    !---------------------------------------------------------------------------

    a11 = a(1, 1, :)
    a21 = a(2, 1, :)
    a22 = a(2, 2, :)
    a31 = a(3, 1, :)
    a32 = a(3, 2, :)
    a33 = a(3, 3, :)
    a41 = a(4, 1, :)
    a42 = a(4, 2, :)
    a43 = a(4, 3, :)
    a44 = a(4, 4, :)
    a51 = a(5, 1, :)
    a52 = a(5, 2, :)
    a53 = a(5, 3, :)
    a54 = a(5, 4, :)
    a55 = a(5, 5, :)

    ! ** a = ldl'.
    ! j = 1

    ai55 = 1.0d0/a11
    a21 = a21*ai55
    a31 = a31*ai55
    a41 = a41*ai55
    a51 = a51*ai55

    ! j = 2

    ai11 = a21*a11
    a22 = a22 - a21*ai11

    a32 = a32 - a31*ai11
    a42 = a42 - a41*ai11
    a52 = a52 - a51*ai11
    ai55 = 1.0d0/a22
    a32 = a32*ai55
    a42 = a42*ai55
    a52 = a52*ai55

    ! j = 3

    ai11 = a31*a11
    ai22 = a32*a22
    ai55 = a31*ai11 + a32*ai22

    a33 = a33 - ai55

    ai55 = 1.0d0/a33
    a43 = a43 - a41*ai11 - a42*ai22
    a53 = a53 - a51*ai11 - a52*ai22
    a43 = a43*ai55
    a53 = a53*ai55

    ! j = 4

    ai11 = a41*a11
    ai22 = a42*a22
    ai33 = a43*a33
    ai55 = a41*ai11 + a42*ai22 + a43*ai33

    a44 = a44 - ai55

    a54 = a54 - a51*ai11 - a52*ai22 - a53*ai33
    a54 = a54/a44

    ! j = 5

    ai11 = a51*a11
    ai22 = a52*a22
    ai33 = a53*a33
    ai44 = a54*a44
    ai55 = a51*ai11 + a52*ai22 + a53*ai33 + a54*ai44

    a55 = a55 - ai55

    ! Column 1 of inverse
    ! Ly = b

    ai21 = -a21
    ai31 = -a31 - a32*ai21
    ai41 = -a41 - a42*ai21 - a43*ai31
    ai51 = -a51 - a52*ai21 - a53*ai31 - a54*ai41

    ! Dz = y

    ai11 = 1.0d0/a11
    ai21 = ai21/a22
    ai31 = ai31/a33
    ai41 = ai41/a44
    ai51 = ai51/a55

    ! l'x = z

    ai41 = ai41 - a54*ai51
    ai31 = ai31 - a43*ai41 - a53*ai51
    ai21 = ai21 - a32*ai31 - a42*ai41 - a52*ai51
    ai11 = ai11 - a21*ai21 - a31*ai31 - a41*ai41 - a51*ai51

    ! Column 2 of inverse
    ! Ly = b

    ai32 = -a32
    ai42 = -a42 - a43*ai32
    ai52 = -a52 - a53*ai32 - a54*ai42

    ! Dz = y

    ai22 = 1.0d0/a22
    ai32 = ai32/a33
    ai42 = ai42/a44
    ai52 = ai52/a55

    ! l'x = z

    ai42 = ai42 - a54*ai52
    ai32 = ai32 - a43*ai42 - a53*ai52
    ai22 = ai22 - a32*ai32 - a42*ai42 - a52*ai52

    ! Column 3 of inverse
    ! Ly = b

    ai43 = -a43
    ai53 = -a53 - a54*ai43

    ! Dz = y

    ai33 = 1.0d0/a33
    ai43 = ai43/a44
    ai53 = ai53/a55

    ! l'x = z

    ai43 = ai43 - a54*ai53
    ai33 = ai33 - a43*ai43 - a53*ai53

    ! Column 4 of inverse
    ! Ly = b

    ai54 = -a54

    ! Dz = y

    ai44 = 1.0d0/a44
    ai54 = ai54/a55

    ! l'x = z

    ai44 = ai44 - a54*ai54

    ! Column 5 of inverse
    ! Dz = y

    ai55 = 1.0d0/a55

    ! Recover Array

    a(1, 1, :) = ai11
    a(2, 1, :) = ai21
    a(3, 1, :) = ai31
    a(4, 1, :) = ai41
    a(5, 1, :) = ai51
    a(2, 2, :) = ai22
    a(3, 2, :) = ai32
    a(4, 2, :) = ai42
    a(5, 2, :) = ai52
    a(3, 3, :) = ai33
    a(4, 3, :) = ai43
    a(5, 3, :) = ai53
    a(4, 4, :) = ai44
    a(5, 4, :) = ai54
    a(5, 5, :) = ai55

    do i = 1, 5
      do j = (i + 1), 5
        a(i, j, :) = a(j, i, :)
      end do
    end do

    return

  end subroutine invert5x5

  !===========================================================================

  subroutine lattice_deform_(qr5x5, vec, vec_lat)

    real(rk), intent(in) :: qr5x5(:, :)
    real(rk), intent(in) :: vec(:)
    real(rk), intent(out) :: vec_lat(:)

    real(rk), allocatable :: qr5x5_b(:,:,:)
    real(rk), allocatable :: vec_b(:,:)
    real(rk), allocatable :: vec_lat_b(:,:)

    integer :: size1, size2

    size1 = size(vec, 1)
    size2 = size(vec_lat, 1)

    allocate(qr5x5_b(5, 5, 1))
    allocate(vec_b(size1, 1))
    allocate(vec_lat_b(size2, 1))

    qr5x5_b(:,:,1) = qr5x5
    vec_b(:, 1) = vec

    call lattice_deform(qr5x5_b, vec_b, vec_lat_b)

    vec_lat = vec_lat_b(:, 1)

  end subroutine lattice_deform_

  subroutine lattice_deform(qr5x5, vec, vec_lat)

    ! Change coordinates from sample to lattice (for 5-vectors).

    ! Rotate{vec}_sam to {vec}_lat via:
    ! {vec_lat} = [qr5x5]' * {vec}

    !----------------------------------------------------------------------

    ! Arguments:
    ! qr5x5: Array of rotation matrices from lattice to sample
    ! vec: Array of 5-vectors in sample basis
    ! vec_lat: Array of 5-vectors in lattice basis
    ! m: Number of elements

    real(rk), intent(in) :: qr5x5(:, :, :)
    real(rk), intent(in) :: vec(:, :)
    real(rk), intent(out) :: vec_lat(:, :)

    integer :: m

    ! Locals:
    ! i, j: Looping indices

    integer :: i
    integer :: j

    m = size(vec, 2)

    !----------------------------------------------------------------------

    ! Rotate {vec}_sm -> {vec}_lat : {v_lat} = [qr5x5]' {v_sm}

    vec_lat = 0.0d0

    do i = 1, 5
      do j = 1, 5
        vec_lat(i, :) = vec_lat(i, :) + qr5x5(j, i, :)*vec(j, :)
      end do
    end do

    return

  end subroutine lattice_deform

  !===========================================================================

  subroutine lattice_spin_(qr3x3, vec, vec_lat)

    real(rk), intent(in) :: qr3x3(3, 3)
    real(rk), intent(in) :: vec(3)
    real(rk), intent(out) :: vec_lat(3)

    real(rk) :: qr3x3_b(3, 3, 1)
    real(rk) :: vec_b(3, 1)
    real(rk) :: vec_lat_b(3, 1)

    qr3x3_b(:,:,1) = qr3x3
    vec_b(:, 1) = vec

    call lattice_spin(qr3x3_b, vec_b, vec_lat_b)

    vec_lat = vec_lat_b(:, 1)

  end subroutine lattice_spin_

  subroutine lattice_spin(qr3x3, w_vec, w_vec_lat)

    ! Convert a skew 3-vector from sample to lattice coordinates.

    ! Rotate {w_vec}_sam to {w_vec}_lat via:
    ! {w_vec_lat} = [qr3x3]' * {w_vec}

    !---------------------------------------------------------------------------

    ! Arguments:
    ! qr3x3: Array of rotation matrices from lattice to sample
    ! w_vec: Array of skew 3-vectors in sample basis
    ! w_vec_lat: Array of skew 3-vectors in lattice basis
    ! m: Number of elements

    real(rk), intent(in) :: qr3x3(:, :, :)
    real(rk), intent(in) :: w_vec(:, :)
    real(rk), intent(out) :: w_vec_lat(:, :)

    ! Locals:
    ! i, j: Looping indices

    integer :: i, j, m

    m = size(qr3x3, 3)
    !---------------------------------------------------------------------------

    w_vec_lat = 0.0d0

    do i = 1, 3
      do j = 1, 3
        w_vec_lat(i, :) = w_vec_lat(i, :) + qr3x3(j, i, :)*w_vec(j, :)
      end do
    end do

    return

  end subroutine lattice_spin

  !===========================================================================

  subroutine mat_vec_skew(mat, vec)

    ! Convert array of skew matrices to array of 3-vectors.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! mat: Array of skew matrices
    ! vec: Array of 3-vectors
    ! m: Number of elements

    real(rk), intent(in) :: mat(:, :, :)
    real(rk), intent(out) :: vec(:, :)

    integer :: m

    m = size(vec, 2)

    !---------------------------------------------------------------------------

    vec = 0.0d0

    vec(1, :) = mat(2, 1, :)
    vec(2, :) = mat(3, 1, :)
    vec(3, :) = mat(3, 2, :)

    return

  end subroutine mat_vec_skew

  !===========================================================================

  subroutine mat_vec_symm_(mat, vec)

    real(rk), intent(in) :: mat(:, :)
    real(rk), intent(out) :: vec(:)

    real(rk) :: mat_b(3, 3, 1)
    real(rk) :: vec_b(5, 1)

    mat_b(:,:,1) = mat
    vec_b(:, 1) = vec

    call mat_vec_symm(mat_b, vec_b)

    vec = vec_b(:, 1)

  end subroutine mat_vec_symm_

  subroutine mat_vec_symm(mat, vec)

    ! Convert array of 3x3 symmetric matrices to array of 5-vectors.

    !----------------------------------------------------------------------

    ! Arguments:
    ! mat: Array of symmetric matrices
    ! vec: Array of 5-vectors
    ! m: Number of elements

    real(rk), intent(in) :: mat(:, :, :)
    real(rk), intent(out) :: vec(:, :)

    ! Locals:
    ! sqr2, sqr32: Local parameters
    integer :: m

    real(rk) :: sqr2
    real(rk) :: sqr32

    m = size(vec, 2)

    !----------------------------------------------------------------------

    sqr2 = dsqrt(2.0d0)
    sqr32 = dsqrt(1.5d0)
    vec = 0.0d0

    vec(1, :) = (mat(1, 1, :) - mat(2, 2, :))/sqr2
    vec(2, :) = mat(3, 3, :)*sqr32
    vec(3, :) = mat(2, 1, :)*sqr2
    vec(4, :) = mat(3, 1, :)*sqr2
    vec(5, :) = mat(3, 2, :)*sqr2

    return

  end subroutine mat_vec_symm

  !===========================================================================

  subroutine mat_x_mat3(a, b, c)

    ! Matrix multiplication for arrays of 3x3 matrices. (c = a*b)

    !---------------------------------------------------------------------------

    ! Arguments:
    ! a, b: Arrays of input matrices
    ! c: Array of output matrices (a*b)

    real(rk), intent(in) :: a(:, :, :)
    real(rk), intent(in) :: b(:, :, :)
    real(rk), intent(out) :: c(:, :, :)

    ! Locals:
    ! i, j, k: Looping indices

    integer :: i
    integer :: j
    integer :: k

    !---------------------------------------------------------------------------

    c = 0.0d0

    do j = 1, 3
      do k = 1, 3
        do i = 1, 3
          c(i, j, :) = c(i, j, :) + a(i, k, :)*b(k, j, :)
        end do
      end do
    end do

    return

  end subroutine mat_x_mat3

  !===========================================================================

  subroutine mat_x_mat5(a, b, c)

    ! Matrix multiplication for arrays of 5x5 matrices.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! a, b: Arrays of 5x5 matrices (input)
    ! c: Array of 5x5 matrices (c=a*b) (output)

    real(rk), intent(in) :: a(:, :, :)
    real(rk), intent(in) :: b(:, :, :)
    real(rk), intent(out) :: c(:, :, :)

    ! Locals:
    ! i, j, k: Looping indices

    integer :: i
    integer :: j
    integer :: k

    !----------------------------------------------------------------------

    c = 0.0d0

    do j = 1, 5
      do k = 1, 5
        do i = 1, 5
          c(i, j, :) = c(i, j, :) + a(i, k, :)*b(k, j, :)
        end do
      end do
    end do

    return

  end subroutine mat_x_mat5

  !===========================================================================

  subroutine mat_x_matt3(a, b, c)

    ! Matrix multiplication by transpose. (3x3)

    !---------------------------------------------------------------------------

    ! Arguments:
    ! a, b: Arrays of 3x3 matrices
    ! c: Array of 3x3 matrices, c=a*b^t

    real(rk), intent(in) :: a(:, :, :)
    real(rk), intent(in) :: b(:, :, :)
    real(rk), intent(out) :: c(:, :, :)

    ! Locals:
    ! i, j, k: Looping indices

    integer :: i
    integer :: j
    integer :: k

    !---------------------------------------------------------------------------

    c = 0.0d0

    do j = 1, 3
      do k = 1, 3
        do i = 1, 3
          c(i, j, :) = c(i, j, :) + a(i, k, :)*b(j, k, :)
        end do
      end do
    end do

    return

  end subroutine mat_x_matt3

  !===========================================================================

  subroutine mat_x_matt5(a, b, c)

    ! Matrix multiplication by transpose for arrays of 5x5 matrices.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! a, b: Arrays of 5x5 matrices (input)
    ! c: Array of 5x5 matrices; c=a*b^t (output)

    real(rk), intent(in) :: a(:, :, :)
    real(rk), intent(in) :: b(:, :, :)
    real(rk), intent(out) :: c(:, :, :)

    ! Locals:
    ! i, j, k: Looping indices

    integer :: i
    integer :: j
    integer :: k

    !---------------------------------------------------------------------------

    c = 0.0d0

    do k = 1, 5
      do j = 1, 5
        do i = 1, 5
          c(i, j, :) = c(i, j, :) + a(i, k, :)*b(j, k, :)
        end do
      end do
    end do

    return

  end subroutine mat_x_matt5

  !===========================================================================

  subroutine mat_x_mats3(a, b, c, m)

    !     Multiply array of 3x3 matrices by a fixed 3x3 matrix, tranposed.

    !----------------------------------------------------------------------

    ! Arguments:
    ! a: Array of matrices to be multiplied
    ! b: Fixed 3x3 matrix
    ! c: Array of 3x3 matrices, c=a*b^t
    ! m: Number of elements

    real(rk), intent(in) :: a(3, 3, m)
    real(rk), intent(in) :: b(3, 3)
    real(rk), intent(out) :: c(3, 3, m)
    integer, intent(in) :: m

    ! Locals:
    ! i, j, k: Looping indices

    integer :: i
    integer :: j
    integer :: k

    !----------------------------------------------------------------------

    c = 0.0d0

    do k = 1, 3
      do j = 1, 3
        do i = 1, 3
          c(i, j, :) = c(i, j, :) + a(i, k, :)*b(j, k)
        end do
      end do
    end do

    return

  end subroutine mat_x_mats3

  !===========================================================================

  subroutine mat_x_vec5(matrix, vector, product)

    ! Multiply array of matrices times array of vectors. (5 dim)

    !----------------------------------------------------------------------

    ! Arguments:
    ! matrix: Array of matrices (5x5)
    ! vector: Array of 5-vectors
    ! product: Array of 5-vectors; product=matrix * vector

    real(rk), intent(in) :: matrix(:, :, :)
    real(rk), intent(in) :: vector(:, :)
    real(rk), intent(out) :: product(:, :)

    ! Locals:
    ! i, j: Looping indices

    integer :: i
    integer :: j

    !----------------------------------------------------------------------

    product = 0.0d0

    do j = 1, 5
      do i = 1, 5
        product(i, :) = product(i, :) + matrix(i, j, :)* &
            & vector(j, :)
      end do
    end do

    return

  end subroutine mat_x_vec5

  !===========================================================================

  subroutine matrix_vec_mult(m, v, mv, dim)

    ! a simple matrix vector multiplication routine that takes advantage of the
    !   compiler's vectorization optimizations. It ends up being faster than
    !   using the built in Fortran matmul method (rc 2017).

    !---------------------------------------------------------------------------

    ! Arguments:
    ! m: Array of matrices
    ! v: Array of vectors
    ! mv: Array of products of m*v
    ! dim:

    real(rk), intent(in) :: m(dim, dim)
    real(rk), intent(in) :: v(dim)
    real(rk), intent(out) :: mv(dim)
    integer, intent(in) :: dim

    ! Locals:

    integer :: i

    !---------------------------------------------------------------------------

    mv = 0.0d0

    do i = 1, dim
      mv(:) = mv(:) + m(:, i)*v(i)
    end do

  end subroutine matrix_vec_mult

  !===========================================================================

  subroutine matt_x_mat3(a, b, c, m)

    ! Matrix multiplication by transpose. (3x3)

    !---------------------------------------------------------------------------

    ! Arguments:
    ! a, b: arrays of 3x3 matrices
    ! c: array of 3x3 matrices, c=a^t*b
    ! m: number of elements

    real(rk), intent(in) :: a(3, 3, m)
    real(rk), intent(in) :: b(3, 3, m)
    real(rk), intent(out) :: c(3, 3, m)
    integer, intent(in) :: m

    ! Locals:
    ! i, j, k: Looping indices

    integer :: i
    integer :: j
    integer :: k

    !----------------------------------------------------------------------

    c = 0.0d0

    do j = 1, 3
      do i = 1, 3
        do k = 1, 3
          c(i, j, :) = c(i, j, :) + a(k, i, :)*b(k, j, :)
        end do
      end do
    end do

    return

  end subroutine matt_x_mat3

  !===========================================================================

  subroutine vec5_norm(vector, norm)

    ! Compute 2-norm for array of 5-vectors.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! norm: Array of vector norms
    ! vector: Array of 5-vectors

    real(rk), intent(in) :: vector(:, :)
    real(rk), intent(out) :: norm(:)

    ! Locals:
    ! i: Looping index

    integer :: i

    !----------------------------------------------------------------------

    norm = 0.0d0

    do i = 1, 5
      norm = norm + vector(i, :)*vector(i, :)
    end do

    norm = dsqrt(norm)

    return

  end subroutine vec5_norm

  !===========================================================================

  subroutine vec3_norm(vector, norm)

    ! Compute 2-norm for array of 3-vectors.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! norm: Array of vector norms
    ! vector: Array of 5-vectors

    real(rk), intent(in) :: vector(:, :)
    real(rk), intent(out) :: norm(:)

    ! Locals:
    ! i: Looping index

    integer :: i

    !----------------------------------------------------------------------

    norm = 0.0d0

    do i = 1, 3
      norm = norm + vector(i, :)*vector(i, :)
    end do

    norm = dsqrt(norm)

    return

  end subroutine vec3_norm

  !===========================================================================

  subroutine rotmat_skew(c, qr3x3)

    ! Construct 3x3 matrix acting on skew 3-vectors which satisfy the matrix
    !   operation:

    ! [w]_sam = [c] * [w]_lat * [c]'

    ! thus, for vectors:

    ! {w}_sam = [qr3x3] * {w}_lat

    !---------------------------------------------------------------------------

    ! Arguments:
    ! c: Array of orientation matrices
    ! qr3x3: Array of 3x3 matrices which acts on skew vectors
    ! m: Number of elements

    real(rk), intent(in)  :: c(:, :, :)
    real(rk), intent(out) :: qr3x3(:, :, :)

    ! Locals:
    ! Cxx: Intermediary calculations

    real(rk), allocatable :: c11(:)
    real(rk), allocatable :: c12(:)
    real(rk), allocatable :: c13(:)
    real(rk), allocatable :: c21(:)
    real(rk), allocatable :: c22(:)
    real(rk), allocatable :: c23(:)
    real(rk), allocatable :: c31(:)
    real(rk), allocatable :: c32(:)
    real(rk), allocatable :: c33(:)

    integer :: m

    m = size(c, 3)

    allocate(c11(m))
    allocate(c12(m))
    allocate(c13(m))
    allocate(c21(m))
    allocate(c22(m))
    allocate(c23(m))
    allocate(c31(m))
    allocate(c32(m))
    allocate(c33(m))

    !---------------------------------------------------------------------------

    c11 = c(1, 1, :)
    c21 = c(2, 1, :)
    c31 = c(3, 1, :)
    c12 = c(1, 2, :)
    c22 = c(2, 2, :)
    c32 = c(3, 2, :)
    c13 = c(1, 3, :)
    c23 = c(2, 3, :)
    c33 = c(3, 3, :)

    qr3x3(1, 1, :) = c22*c11 - c21*c12
    qr3x3(1, 2, :) = c23*c11 - c21*c13
    qr3x3(1, 3, :) = c23*c12 - c22*c13
    qr3x3(2, 1, :) = c32*c11 - c31*c12
    qr3x3(2, 2, :) = c33*c11 - c31*c13
    qr3x3(2, 3, :) = c33*c12 - c32*c13
    qr3x3(3, 1, :) = c32*c21 - c31*c22
    qr3x3(3, 2, :) = c33*c21 - c31*c23
    qr3x3(3, 3, :) = c33*c22 - c32*c23

    return

  end subroutine rotmat_skew

  subroutine rotmat_skew_(c, qr3x3)

    real(rk), intent(in)  :: c(:, :)
    real(rk), intent(out) :: qr3x3(:, :)

    qr3x3(1, 1) = c(2, 2)*c(1, 1) - c(2, 1)*c(1, 2)
    qr3x3(1, 2) = c(2, 3)*c(1, 1) - c(2, 1)*c(1, 3)
    qr3x3(1, 3) = c(2, 3)*c(1, 2) - c(2, 2)*c(1, 3)
    qr3x3(2, 1) = c(3, 2)*c(1, 1) - c(3, 1)*c(1, 2)
    qr3x3(2, 2) = c(3, 3)*c(1, 1) - c(3, 1)*c(1, 3)
    qr3x3(2, 3) = c(3, 3)*c(1, 2) - c(3, 2)*c(1, 3)
    qr3x3(3, 1) = c(3, 2)*c(2, 1) - c(3, 1)*c(2, 2)
    qr3x3(3, 2) = c(3, 3)*c(2, 1) - c(3, 1)*c(2, 3)
    qr3x3(3, 3) = c(3, 3)*c(2, 2) - c(3, 2)*c(2, 3)

    return

  end subroutine rotmat_skew_

  subroutine rotmat_skew2(c, qr3x3, m)

    ! Construct 3x3 matrix acting on (properly ordered) rotation vector which satisfy the matrix
    !   operation:

    ! [w]_sam = [c] * [w]_lat * [c]'

    ! thus, for vectors:

    ! {w}_sam = [qr3x3] * {w}_lat

    !---------------------------------------------------------------------------

    ! Arguments:
    ! c: Array of orientation matrices
    ! qr3x3: Array of 3x3 matrices
    ! m: Number of elements

    real(rk), intent(in)  :: c(3, 3, m)
    real(rk), intent(out) :: qr3x3(3, 3, m)
    integer, intent(in) :: m

    ! Locals:
    ! Cxx: Intermediary calculations

    real(rk) :: c11(m)
    real(rk) :: c12(m)
    real(rk) :: c13(m)
    real(rk) :: c21(m)
    real(rk) :: c22(m)
    real(rk) :: c23(m)
    real(rk) :: c31(m)
    real(rk) :: c32(m)
    real(rk) :: c33(m)

    !---------------------------------------------------------------------------

    ! permutation of indices 1 and 3. -1 on the 2 index
    c11 = c(3, 3, :)
    c21 = -c(2, 3, :)
    c31 = c(1, 3, :)
    c12 = -c(3, 2, :)
    c22 = c(2, 2, :)
    c32 = -c(1, 2, :)
    c13 = c(3, 1, :)
    c23 = -c(2, 1, :)
    c33 = c(1, 1, :)

    qr3x3(1, 1, :) = c22*c11 - c21*c12
    qr3x3(1, 2, :) = c23*c11 - c21*c13
    qr3x3(1, 3, :) = c23*c12 - c22*c13
    qr3x3(2, 1, :) = c32*c11 - c31*c12
    qr3x3(2, 2, :) = c33*c11 - c31*c13
    qr3x3(2, 3, :) = c33*c12 - c32*c13
    qr3x3(3, 1, :) = c32*c21 - c31*c22
    qr3x3(3, 2, :) = c33*c21 - c31*c23
    qr3x3(3, 3, :) = c33*c22 - c32*c23

    return

  end subroutine rotmat_skew2

  !===========================================================================

  subroutine rotmat_symm(c, qr5x5, m)

    ! Convert rotation matrix c to an operator on 5-vectors which satisfy the
    !   matrix operation:

    ! [a]_sam = [c] * [a]_lat * [c]'

    ! thus, for vectors:

    ! {a}_sam = [qr5x5] * {a}_lat

    ! with: {a}={()/sqr2,sqr32*(),sqr2*(),sqr2*(),sqr2*()}

    !---------------------------------------------------------------------------

    ! Arguments:
    ! c: Array of rotation matrices
    ! qr5x5: Array of 5x5 rotation matrices
    ! m: Number of elements

    real(rk), intent(in)  :: c(3, 3, m)
    real(rk), intent(out) :: qr5x5(5, 5, m)
    integer, intent(in) :: m

    ! Locals:
    ! sqr3: Local parameter
    ! Cxx: Intermediary calculations

    real(rk) :: sqr3
    real(rk) :: c11(m)
    real(rk) :: c12(m)
    real(rk) :: c13(m)
    real(rk) :: c21(m)
    real(rk) :: c22(m)
    real(rk) :: c23(m)
    real(rk) :: c31(m)
    real(rk) :: c32(m)
    real(rk) :: c33(m)

    !----------------------------------------------------------------------

    sqr3 = dsqrt(3.0d0)

    c11 = c(1, 1, :)
    c21 = c(2, 1, :)
    c31 = c(3, 1, :)
    c12 = c(1, 2, :)
    c22 = c(2, 2, :)
    c32 = c(3, 2, :)
    c13 = c(1, 3, :)
    c23 = c(2, 3, :)
    c33 = c(3, 3, :)

    qr5x5(1, 1, :) = 0.5d0*(c11*c11 - c12*c12 - c21*c21 + c22* &
        & c22)
    qr5x5(1, 2, :) = sqr3/2.0d0*(c13*c13 - c23*c23)
    qr5x5(1, 3, :) = c11*c12 - c21*c22
    qr5x5(1, 4, :) = c11*c13 - c21*c23
    qr5x5(1, 5, :) = c12*c13 - c22*c23
    qr5x5(2, 1, :) = sqr3/2.0d0*(c31*c31 - c32*c32)
    qr5x5(2, 2, :) = 1.5d0*c33*c33 - 0.5d0
    qr5x5(2, 3, :) = sqr3*c31*c32
    qr5x5(2, 4, :) = sqr3*c31*c33
    qr5x5(2, 5, :) = sqr3*c32*c33
    qr5x5(3, 1, :) = c11*c21 - c12*c22
    qr5x5(3, 2, :) = sqr3*c13*c23
    qr5x5(3, 3, :) = c11*c22 + c12*c21
    qr5x5(3, 4, :) = c11*c23 + c13*c21
    qr5x5(3, 5, :) = c12*c23 + c13*c22
    qr5x5(4, 1, :) = c11*c31 - c12*c32
    qr5x5(4, 2, :) = sqr3*c13*c33
    qr5x5(4, 3, :) = c11*c32 + c12*c31
    qr5x5(4, 4, :) = c11*c33 + c13*c31
    qr5x5(4, 5, :) = c12*c33 + c13*c32
    qr5x5(5, 1, :) = c21*c31 - c22*c32
    qr5x5(5, 2, :) = sqr3*c23*c33
    qr5x5(5, 3, :) = c21*c32 + c22*c31
    qr5x5(5, 4, :) = c21*c33 + c23*c31
    qr5x5(5, 5, :) = c22*c33 + c23*c32

    return

  end subroutine rotmat_symm

  !===========================================================================

  subroutine solve_lin_sys_3(mat, vec, sol)

    ! Solve a linear system of three equations. Used for triaxial loading.

    !---------------------------------------------------------------------------

    ! Arguments:

    real(rk), intent(in) :: mat(3, 3)
    real(rk), intent(in) :: vec(3)
    real(rk), intent(out) :: sol(3)

    ! Locals:

    real(rk) :: inv(3, 3)
    real(rk) :: a, b, c, d, e, f, g, h, k
    real(rk) :: det, matnorm, invnorm, cond

    !---------------------------------------------------------------------------

    a = mat(1, 1)
    b = mat(1, 2)
    c = mat(1, 3)
    d = mat(2, 1)
    e = mat(2, 2)
    f = mat(2, 3)
    g = mat(3, 1)
    h = mat(3, 2)
    k = mat(3, 3)

    det = a*(e*k - f*h) - b*(d*k - f*g) + c*(d*h - e*g)

    inv(1, 1) = (e*k - f*h)/det
    inv(1, 2) = (c*h - b*k)/det
    inv(1, 3) = (b*f - c*e)/det
    inv(2, 1) = (f*g - d*k)/det
    inv(2, 2) = (a*k - c*g)/det
    inv(2, 3) = (c*d - a*f)/det
    inv(3, 1) = (d*h - e*g)/det
    inv(3, 2) = (b*g - a*h)/det
    inv(3, 3) = (a*e - b*d)/det

    sol = matmul(inv, vec)

    ! Find conditioning number

    matnorm = max(mat(1, 1) + mat(1, 2) + mat(1, 3), &
         & mat(2, 1) + mat(2, 2) + mat(2, 3), &
         & mat(3, 1) + mat(3, 2) + mat(3, 3))
    invnorm = max(inv(1, 1) + inv(1, 2) + inv(1, 3), &
         & inv(2, 1) + inv(2, 2) + inv(2, 3), &
         & inv(3, 1) + inv(3, 2) + inv(3, 3))
    cond = matnorm*invnorm

    if (cond .gt. 1.0d3) then
      call par_quit('Error  :     > Matrix is poorly conditioned.')
    end if

    return

  end subroutine solve_lin_sys_3

  !===========================================================================

  subroutine sparse_matvec_ebe(res, sol, temp1, temp2, gstif, bcs_vel_defined, kdim, &
      & dof_sub, dof_sup, elt_sub, elt_sup, dtrace, elt_dofs)

    !  Matrix vector multiply, element by element.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! res: 
    ! sol: Values on the global nodal array
    ! temp1: Values from sol array gather in elemental array using connectivity
    ! temp2:
    ! gstif: Global stiffness matrix (collection of elemental stiffness matrix)
    ! bcs_vel_defined: Boundary condition (true/false) array
    ! kdim: Dimension of the elemental stifness matrix (kdim*kdim)
    ! dof_sub: min dof for this process
    ! dof_sup: max dof for this process
    ! elt_sub: min elt for this process
    ! elt_sup: max elt for this process
    ! dtrace: Communication information for gather/scatter
    ! elt_dofs: Connectivity

    real(rk) :: res(dof_sub:dof_sup)
    real(rk) :: sol(dof_sub:dof_sup)
    real(rk) :: temp1(kdim, elt_sub:elt_sup)
    real(rk) :: temp2(kdim, elt_sub:elt_sup)
    real(rk) :: gstif(kdim, kdim, elt_sub:elt_sup)
    logical :: bcs_vel_defined(dof_sub:dof_sup)
    integer :: kdim
    integer :: dof_sub
    integer :: dof_sup
    integer :: elt_sub
    integer :: elt_sup
    type(trace) :: dtrace
    integer :: elt_dofs(ndim, elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    call part_gather(temp1, sol, elt_dofs, dtrace)

    call gen_matrix_vector_mult(gstif, temp1, temp2)

    call part_scatter(res, temp2, elt_dofs, dtrace)

    where (bcs_vel_defined)
      res = 0.0d0
    end where

    return

  end subroutine sparse_matvec_ebe

  !===========================================================================

  subroutine velgrad_sympart(velgrad, d, dkk)

    ! Compute symmetric part of an array of vel gradients.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! d: Deviatoric part of vel gradients
    ! dkk: Trace of vel gradients
    ! velgrad: Array of vel gradients

    real(rk), intent(in) :: velgrad(:, :, :)
    real(rk), intent(out) :: d(:, :, :)
    real(rk), intent(out) :: dkk(:)

    !---------------------------------------------------------------------------

    d(1, 1, :) = velgrad(1, 1, :)
    d(2, 2, :) = velgrad(2, 2, :)
    d(3, 3, :) = velgrad(3, 3, :)
    d(2, 1, :) = 0.5d0*(velgrad(2, 1, :) + velgrad(1, 2, :))
    d(3, 1, :) = 0.5d0*(velgrad(3, 1, :) + velgrad(1, 3, :))
    d(3, 2, :) = 0.5d0*(velgrad(3, 2, :) + velgrad(2, 3, :))
    d(1, 2, :) = d(2, 1, :)
    d(1, 3, :) = d(3, 1, :)
    d(2, 3, :) = d(3, 2, :)

    dkk = d(1, 1, :) + d(2, 2, :) + d(3, 3, :)

    d(1, 1, :) = d(1, 1, :) - dkk/3.0d0
    d(2, 2, :) = d(2, 2, :) - dkk/3.0d0
    d(3, 3, :) = d(3, 3, :) - dkk/3.0d0

    return

  end subroutine velgrad_sympart

  !===========================================================================

  subroutine velgrad_skewpart(velgrad, w)

    ! Compute skew part of an array of vel gradients.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! w: Array of skew matrices (output)
    ! velgrad: Array of vel gradients
    ! m: Number of elements

    real(rk), intent(in) :: velgrad(:, :, :)
    real(rk), intent(out) :: w(:, :, :)

    !---------------------------------------------------------------------------

    w(:, :, :) = 0.0d0

    w(1, 2, :) = 0.5d0*(velgrad(1, 2, :) - velgrad(2, 1, :))
    w(1, 3, :) = 0.5d0*(velgrad(1, 3, :) - velgrad(3, 1, :))
    w(2, 3, :) = 0.5d0*(velgrad(2, 3, :) - velgrad(3, 2, :))
    w(2, 1, :) = -w(1, 2, :)
    w(3, 1, :) = -w(1, 3, :)
    w(3, 2, :) = -w(2, 3, :)

    return

  end subroutine velgrad_skewpart

  !===========================================================================

  subroutine strain_equiv_3x3(e, equiv)

    ! Find the equivalent value of a strain tensor,
    !   equiv = sqrt(2/3 e_(ij) e_(ij)).

    !---------------------------------------------------------------------------

    ! Arguments:
    ! e: Array of 3x3 matrices

    real(rk), intent(in) :: e(3, 3, elt_sub:elt_sup)
    real(rk), intent(out) :: equiv(elt_sub:elt_sup)

    ! Locals:
    ! i: Looping index
    ! edev: Deviatoric portion of e
    ! evol: Volumetric portion of e

    integer :: i
    real(rk) :: evol(elt_sub:elt_sup)
    real(rk) :: edev(3, 3, elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    equiv = 0.0d0
    evol = 0.0d0
    edev = 0.0d0

    evol(:) = (1.0d0/3.0d0)*(e(1, 1, :) + e(2, 2, :) + e(3, 3, :))
    edev = e
    edev(1, 1, :) = edev(1, 1, :) - evol(:)
    edev(2, 2, :) = edev(2, 2, :) - evol(:)
    edev(3, 3, :) = edev(3, 3, :) - evol(:)

    do i = elt_sub, elt_sup
      equiv(i) = dsqrt((2.0d0/3.0d0)*(edev(1, 1, i)*edev(1, 1, i) + &
          & edev(1, 2, i)*edev(1, 2, i) + edev(1, 3, i)*edev(1, 3, i) + &
          & edev(2, 1, i)*edev(2, 1, i) + edev(2, 2, i)*edev(2, 2, i) + &
          & edev(2, 3, i)*edev(2, 3, i) + edev(3, 1, i)*edev(3, 1, i) + &
          & edev(3, 2, i)*edev(3, 2, i) + edev(3, 3, i)*edev(3, 3, i)))
    end do

  end subroutine strain_equiv_3x3

  !> Compute the effective deformation rate from the deformation rate tensor
  subroutine defrate_eq(d, eq)

    real(rk), intent(out) :: eq(elt_sub:elt_sup)
    real(rk), intent(in) :: d(3, 3, elt_sub:elt_sup)

    real(rk), parameter :: twothirds = 2.0d0/3.0d0
    integer :: i

    !---------------------------------------------------------------------------

    do i = elt_sub, elt_sup
      eq(i) = dsqrt(twothirds*(d(1, 1, i)*d(1, 1, i) + d(2, 2, i)* &
          & d(2, 2, i) + d(3, 3, i)*d(3, 3, i) + 2.0d0*(d(1, 2, i)* &
          & d(1, 2, i) + d(1, 3, i)*d(1, 3, i) + d(2, 3, i)*d(2, 3, i))))
    end do

  end subroutine defrate_eq

  subroutine strain_equiv(e, equiv)

    ! Find the equivalent value of a strain tensor,
    !   equiv = sqrt(2/3 e_(ij) e_(ij)).

    !---------------------------------------------------------------------------

    ! Arguments:
    ! e: Array of 3x3 matrices

    real(rk), intent(in) :: e(:, :)
    real(rk), intent(out) :: equiv(:)

    real(rk), parameter :: twothirds = 2.0d0/3.0d0

    equiv = dsqrt(twothirds*sum(e*e, 1))

  end subroutine strain_equiv

  subroutine vector_tiny_zero (vector_in, vector_out)

    real(rk), intent(in) :: vector_in (:)
    real(rk), intent(out) :: vector_out (:)

    vector_out = vector_in

    where (vector_out .le. vtiny) vector_out = 0.0d0

  end subroutine vector_tiny_zero

  !===========================================================================

  subroutine stress_equiv_3x3(s, equiv)

    ! Find the equivalent value of a strain tensor,
    !   equiv = sqrt(3/2 s_(ij) s_(ij)).

    !---------------------------------------------------------------------------

    ! Arugments:
    ! s3x3: Stress tensor in 3x3 matrix form
    ! equiv: Equivalent (von Mises) stress (scalar)

    real(rk), intent(in) :: s(3, 3, elt_sub:elt_sup)
    real(rk), intent(out) :: equiv(elt_sub:elt_sup)

    ! Locals:
    ! i: Looping index

    integer :: i
    real(rk) :: svol(elt_sub:elt_sup)
    real(rk) :: sdev(3, 3, elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    equiv = 0.0d0
    svol = 0.0d0
    sdev = 0.0d0

    svol(:) = (1.0d0/3.0d0)*(s(1, 1, :) + s(2, 2, :) + s(3, 3, :))
    sdev = s
    sdev(1, 1, :) = sdev(1, 1, :) - svol(:)
    sdev(2, 2, :) = sdev(2, 2, :) - svol(:)
    sdev(3, 3, :) = sdev(3, 3, :) - svol(:)

    do i = elt_sub, elt_sup
      equiv(i) = dsqrt((3.0d0/2.0d0)*(sdev(1, 1, i)*sdev(1, 1, i) + &
          & sdev(1, 2, i)*sdev(1, 2, i) + sdev(1, 3, i)*sdev(1, 3, i) + &
          & sdev(2, 1, i)*sdev(2, 1, i) + sdev(2, 2, i)*sdev(2, 2, i) + &
          & sdev(2, 3, i)*sdev(2, 3, i) + sdev(3, 1, i)*sdev(3, 1, i) + &
          & sdev(3, 2, i)*sdev(3, 2, i) + sdev(3, 3, i)*sdev(3, 3, i)))
    end do

  end subroutine stress_equiv_3x3

  !===========================================================================

  subroutine tensor3dcompose(mat, dev, skw, sph)

    ! Form  matrix from deviatoric, skew, or spherical parts. If no parts are
    !   passed, the matix is zeroed.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! mat: Resulting array of 3x3 matrices
    ! dev: Arry of 5-vec representing symmetric, traceless portion
    ! skw: Array of 3-vec representing axial vectors of the skew pportion
    ! sph: Array of scalars representing one third of the trace

    real(rk), intent(out) :: mat(:, :, :)
    real(rk), intent(in), optional :: dev(:, :)
    real(rk), intent(in), optional :: skw(:, :)
    real(rk), intent(in), optional :: sph(:)

    !---------------------------------------------------------------------------

    mat = 0.0d0

    if (present(dev)) then
      mat(1, 1, :) = -dev(1, :)*sq2_i - dev(2, :)*sq6_i
      mat(2, 2, :) = dev(1, :)*sq2_i - dev(2, :)*sq6_i
      mat(3, 3, :) = dev(2, :)*twosq6_i

      mat(2, 3, :) = dev(3, :)*sq2_i
      mat(3, 2, :) = mat(2, 3, :)

      mat(1, 3, :) = dev(4, :)*sq2_i
      mat(3, 1, :) = mat(1, 3, :)

      mat(1, 2, :) = dev(5, :)*sq2_i
      mat(2, 1, :) = mat(1, 2, :)
    end if

    if (present(skw)) then
      mat(2, 3, :) = mat(2, 3, :) - skw(1, :)
      mat(3, 2, :) = mat(3, 2, :) + skw(1, :)

      mat(1, 3, :) = mat(1, 3, :) + skw(2, :)
      mat(3, 1, :) = mat(3, 1, :) - skw(2, :)

      mat(1, 2, :) = mat(1, 2, :) - skw(3, :)
      mat(2, 1, :) = mat(2, 1, :) + skw(3, :)
    end if

    if (present(sph)) then
      mat(1, 1, :) = mat(1, 1, :) + sph
      mat(2, 2, :) = mat(2, 2, :) + sph
      mat(3, 3, :) = mat(3, 3, :) + sph
    end if

  end subroutine tensor3dcompose

  !===========================================================================

  subroutine vec5_vec6(vec, mat)

    ! Convert 5-vector to 6-vector of symmetric matrix.

    ! Note: Returns the upper triangle in the order of: 11 22 33 23 13 12

    !---------------------------------------------------------------------------

    ! Arguments:
    ! vec: 5-vector of dev. part of 2nd-order tensor (inner product preserved)
    ! mat: 6-vector of upper triangle of symmetric tensor (no ip preservation)

    real(rk), intent(in) :: vec(5)
    real(rk), intent(out) :: mat(6)

    ! Locals:
    ! sqr2, sqr23: Local parameters

    real(rk) :: sqr2
    real(rk) :: sqr23

    !---------------------------------------------------------------------------

    sqr2 = dsqrt(2.0d0)
    sqr23 = dsqrt(2.0d0/3.0d0)

    ! Diagonal: 11 22 33
    mat(1) = 0.5d0*(sqr2*vec(1) - sqr23*vec(2))
    mat(2) = -0.5d0*(sqr2*vec(1) + sqr23*vec(2))
    mat(3) = vec(2)*sqr23

    ! Off-diagonal: 23 13 12
    mat(4) = vec(5)/sqr2
    mat(5) = vec(4)/sqr2
    mat(6) = vec(3)/sqr2

    return

  end subroutine vec5_vec6

  !===========================================================================

  subroutine vec6_crystosam(vec6_crys, ori, vec6_sam)

    ! Convert the 6 vector of in the crystal basis to the sample basis.

    !---------------------------------------------------------------------------

    ! Arugments:
    ! vec6_crys: 6-vector in the crystal basis
    ! ori: 3x3 orientation matrix in the crystal to sample convention
    ! vec6_sam: 6-vector in the sample basis
    real(rk), intent(in) :: vec6_crys(6, elt_sub:elt_sup)
    real(rk), intent(in) :: ori(3, 3, elt_sub:elt_sup)
    real(rk), intent(out) :: vec6_sam(6, elt_sub:elt_sup)
    ! Locals:
    ! i, j, k, ii: Looping indices
    ! mat33_crys: 3x3 matrix in crystal frame
    ! temp: Matrix of strain*r'
    ! mat33_sam: 3x3 matrix in sample frame
    integer :: i, j, k, ii
    real(rk) :: mat33_crys(3, 3)
    real(rk) :: temp(3, 3)
    real(rk) :: mat33_sam(3, 3)

    !---------------------------------------------------------------------------

    vec6_sam = 0.0d0

    do ii = elt_sub, elt_sup

      mat33_crys = 0.0d0
      mat33_sam = 0.0d0
      temp = 0.0d0

      ! Put 6-vector into 3x3 matrix
      mat33_crys(1, 1) = vec6_crys(1, ii)
      mat33_crys(1, 2) = vec6_crys(2, ii)
      mat33_crys(1, 3) = vec6_crys(3, ii)
      mat33_crys(2, 1) = vec6_crys(2, ii)
      mat33_crys(2, 2) = vec6_crys(4, ii)
      mat33_crys(2, 3) = vec6_crys(5, ii)
      mat33_crys(3, 1) = vec6_crys(3, ii)
      mat33_crys(3, 2) = vec6_crys(5, ii)
      mat33_crys(3, 3) = vec6_crys(6, ii)

      ! Perform strain*r'
      do j = 1, 3
        do k = 1, 3
          do i = 1, 3
            temp(i, j) = temp(i, j) + mat33_crys(i, k)*ori(j, k, ii)
          end do
        end do
      end do

      ! Perform r*(strain*r')
      do j = 1, 3
        do k = 1, 3
          do i = 1, 3
            mat33_sam(i, j) = mat33_sam(i, j) + ori(i, k, ii)*temp(k, j)
          end do
        end do
      end do

      ! Put into 6-vector of order 11 22 33 23 13 12

      vec6_sam(1, ii) = mat33_sam(1, 1)
      vec6_sam(2, ii) = mat33_sam(2, 2)
      vec6_sam(3, ii) = mat33_sam(3, 3)
      vec6_sam(4, ii) = mat33_sam(2, 3)
      vec6_sam(5, ii) = mat33_sam(1, 3)
      vec6_sam(6, ii) = mat33_sam(1, 2)

    end do

  end subroutine vec6_crystosam

  !===========================================================================

  subroutine vec_d_vec5(d, b, c, m)

    ! Computes c = d*b where c and b are 5-vectors and d is a diagonal 5x5
    !   matrix represented in compact form

    !----------------------------------------------------------------------

    ! Arguments:
    ! d: a fixed diagonal matrix (in 5-vector form)
    ! b: Input matrix (in 5-vector form)
    ! c: Output matrix (d*b) (in 5-vector form)
    ! m: Number of elements

    real(rk), intent(in) :: d(5)
    real(rk), intent(in) :: b(5, m)
    real(rk), intent(out) :: c(5, m)
    integer, intent(in) :: m

    ! Locals:
    ! k: Looping index

    integer :: k

    !----------------------------------------------------------------------

    do k = 1, 5
      c(k, :) = d(k)*b(k, :)
    end do

    return

  end subroutine vec_d_vec5

  !===========================================================================

  subroutine vec_mat_skew(vec, mat, m)

    ! Convert array of 3-vectors to array of skew matrices.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! vec: Array of 3-vectors
    ! mat: Array of skew matrices
    ! m: Number of elements

    real(rk), intent(in) :: vec(3, m)
    real(rk), intent(out) :: mat(3, 3, m)
    integer, intent(in) :: m

    !---------------------------------------------------------------------------

    mat(1, 1, :) = 0.0d0
    mat(2, 2, :) = 0.0d0
    mat(3, 3, :) = 0.0d0

    mat(2, 1, :) = vec(1, :)
    mat(3, 1, :) = vec(2, :)
    mat(3, 2, :) = vec(3, :)

    mat(1, 2, :) = -vec(1, :)
    mat(1, 3, :) = -vec(2, :)
    mat(2, 3, :) = -vec(3, :)

    return

  end subroutine vec_mat_skew

  !===========================================================================

  subroutine vec_mat_symm(vec, mat)

    ! Convert 5-vector to symmetric matrix: {5} --> [3x3]sym

    !---------------------------------------------------------------------------

    ! Arguments:
    ! vec: Array of 5-vectors
    ! mat: Array of 3x3 symmetric matrices
    ! m: number of elements

    real(rk), intent(in) :: vec(:, :)
    real(rk), intent(out) :: mat(:, :, :)
    integer :: m

    ! Locals:
    ! sqr2, sqr23: Local parameters

    real(rk) :: sqr2
    real(rk) :: sqr23

    !----------------------------------------------------------------------

    m = size(vec, 2)

    sqr2 = dsqrt(2.0d0)
    sqr23 = dsqrt(2.0d0/3.0d0)

    mat(1, 1, :) = 0.5d0*(sqr2*vec(1, :) - sqr23*vec(2, :))
    mat(2, 2, :) = -0.5d0*(sqr2*vec(1, :) + sqr23*vec(2, :))
    mat(3, 3, :) = vec(2, :)*sqr23
    mat(2, 1, :) = vec(3, :)/sqr2
    mat(3, 1, :) = vec(4, :)/sqr2
    mat(3, 2, :) = vec(5, :)/sqr2

    mat(1, 2, :) = mat(2, 1, :)
    mat(1, 3, :) = mat(3, 1, :)
    mat(2, 3, :) = mat(3, 2, :)

    return

  end subroutine vec_mat_symm

  !===========================================================================

  subroutine vec6_mat_symm(vec, mat, m)

    ! Convert 6-vector to symmetric matrix: {6} --> [3x3]sym

    !---------------------------------------------------------------------------

    ! Arguments:
    ! vec: Array of 6-vectors
    ! mat: Array of 3x3 symmetric matrices
    ! m: number of elements

    real(rk), intent(in) :: vec(6, m)
    real(rk), intent(out) :: mat(3, 3, m)
    integer, intent(in) :: m

    ! Locals:

    !----------------------------------------------------------------------

    mat(1, 1, :) = vec(1, :)
    mat(1, 2, :) = vec(2, :)
    mat(1, 3, :) = vec(3, :)
    mat(2, 1, :) = vec(2, :)
    mat(2, 2, :) = vec(4, :)
    mat(2, 3, :) = vec(5, :)
    mat(3, 1, :) = vec(3, :)
    mat(3, 2, :) = vec(5, :)
    mat(3, 3, :) = vec(6, :)

    return

  end subroutine vec6_mat_symm

  subroutine mat_vec6 (mat, vec6)

    real(rk), intent(in) :: mat (:, :, :)
    real(rk), intent(out) :: vec6 (:, :)

    vec6(1, :) = mat(1, 1, :)
    vec6(2, :) = mat(1, 2, :)
    vec6(3, :) = mat(1, 3, :)
    vec6(4, :) = mat(2, 2, :)
    vec6(5, :) = mat(2, 3, :)
    vec6(6, :) = mat(3, 3, :)

  end subroutine mat_vec6

  !===========================================================================

end module matrix_operations_mod
