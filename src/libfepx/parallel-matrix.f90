! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE parallel_matrix_mod
  !
  !  Parallel matrix ops.
  !
  USE IntrinsicTypesModule, RK=>REAL_KIND
  !
  IMPLICIT NONE
  !
CONTAINS
  !
  !  gen_matrix_vector_mult : matrix times vector
  !  gen_matrix_mult        : matrix time matrix
  !  sparse_matvec_ebe      : matrix times vector, element by element
  !
  ! **********************************************************************
  !
  SUBROUTINE sparse_matvec_ebe(res, sol, temp1, temp2, gstif, bcs, &
       &  nnpe, nsub, nsup, esub, esup, dtrace, np)
    !
    !  Matrix vector multiply, element by element.
    !
    !----------------------------------------------------------------------
    !
    USE gather_scatter
    !
    IMPLICIT NONE
    !
    !  Arguments.
    !
    TYPE(trace) dtrace
    !
    INTEGER nnpe, nsub, nsup, esub, esup
    INTEGER   np(0:(nnpe-1),esub:esup)
    !    
    LOGICAL bcs(nsub:nsup)
    !
    REAL(RK) res(nsub:nsup), sol(nsub:nsup)
    REAL(RK) temp1(0:(nnpe - 1), esub:esup)
    REAL(RK) temp2(0:(nnpe - 1), esub:esup)
    REAL(RK) gstif(0:(nnpe - 1), 0:(nnpe - 1), esub:esup)
    !
    !  Locals:
    !
    INTEGER ier, i, j, idummy
    !
    !----------------------------------------------------------------------
    !
    call part_gather(temp1, sol, np, dtrace)

    temp2 = 0.0_RK
      
    call gen_matrix_vector_mult(temp2, gstif, temp1, &
         &  idummy, idummy, idummy, idummy, ier)
      
    res = 0.0_RK

    call part_scatter(res, temp2, np, .false., dtrace)
 
    where (bcs) res = 0.0_RK
      
    RETURN
  END SUBROUTINE sparse_matvec_ebe
  !
  ! *********************************************************************
  !
  SUBROUTINE gen_matrix_vector_mult(y, a, x, i1, i2, i3, i4, ier)
    !
    !  Array of matrices times array of vectors: y(i) = A(i)*x(i)
    !
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !  Arguments:
    !
    REAL(RK) y(:,:), x(:,:), a(:,:,:)
    INTEGER i1, i2, i3, i4, ier
    !
    !  Locals:
    !
    INTEGER i, j, k
    INTEGER ad1, ad2, ad3
    !
    !----------------------------------------------------------------------
    !
    !     Note that the arrays x,y,a are all assumed shape arrays and
    !     that they are therefore dimensioned 1:p, 1:q, etc., even though
    !     they may be dimensioned differently in the calling routine.
    !
    ad1 = ubound(a, 1)
    ad2 = ubound(a, 2)
    ad3 = ubound(a, 3)
    !
    do k=1, ad3
       do j=1,ad2
!          do i=1,ad1
             y(:,k) = y(:,k) + a(:,j,k)*x(j,k)
!          enddo
       enddo
    enddo
    ier = 0
    !
    RETURN
  END SUBROUTINE gen_matrix_vector_mult
  !
  ! **********************************************************************
  !
  SUBROUTINE gen_matrix_mult(c, a, b, i1, i2, ier)
    !
    !     Accumulate Array of matrices times array of matrices:
    !     *     c(i) = a(i) * b(i)
    !
    !     Note: c is not zeroed in this subroutine, so the
    !     *     values are accumulated; in all the calling
    !     *     routines, however, the c array is zeroed before
    !     *     calling this routine.
    !
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !     Arguments:
    !
    !     c     : array of result matrices (c=a*b)
    !     a,b   : array of input matrices
    !     i1,i2 : not used
    !     ier   : return status, 0=okay
    !     
    REAL(RK) c(:,:,:), a(:,:,:), b(:,:,:)
    INTEGER i1, i2, ier
    !
    !     Locals:
    !
    INTEGER i, j, k, m
    INTEGER lda, ldb, ldc, lta, ltb, n
    !
    !----------------------------------------------------------------------
    !
    lda = ubound(a, 1)
    ldb = ubound(b, 1)
    ldc = ubound(c, 1)
    lta = ubound(a, 2)
    ltb = ubound(b, 2)
    n   = ubound(a, 3)
    !
    do m=1, n
       do j=1,ltb
          do k=1,lta
!             do i=1,ldc
                c(:,j,m) = c(:,j,m) + a(:,k,m)*b(k,j,m)
!             enddo
          enddo
       enddo
    enddo
    !
    ier = 0
    !
    RETURN
  END SUBROUTINE gen_matrix_mult
  !
END MODULE parallel_matrix_mod
!
