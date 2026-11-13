! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Mandel 6-vector algebra of symmetric second-order tensors, skew 3-vectors, the 6x6
!> operators acting on Mandel vectors (rotation [Q6(R)], spin commutator [W^]), and deviatoric
!> 5-vectors with the 5x5 form of 6x6 operators.
!>
!> Theory : dev/conventions.md §1.1 (Mandel 6-vectors, order 11 22 33 23 13 12, sqrt(2) on the
!>          shear components), §1.2 (tensor-component 6-vectors for I/O), §1.4 (skew 3-vectors,
!>          order w23 w13 w12); dev/theory/recoupled_formulation.md Eq. (R4) ([Q6(R)]) and
!>          Eq. (R8) ([W^]). The Mandel (orthonormal-basis) representation of second- and
!>          fourth-order tensors: Mehrabadi, M.M. and Cowin, S.C. (1990), "Eigentensors of
!>          linear anisotropic elastic materials", Q. J. Mech. Appl. Math. 43(1):15-41.
!>          Deviatoric 5-vectors: dev/theory/addendum_5vectors.md, after Kasemer, M. and
!>          Dawson, P.R. (2025), arXiv:2502.19531 (KD25), Eqs. 7-9 and Appendix A Eqs. 28-31.
!> Source : equivalent measures as v2 matrix_operations_mod.f90::stress_equiv_3x3 (von Mises
!>          stress, sqrt(3/2 s':s')) and ::strain_equiv_3x3 / ::defrate_eq / ::strain_equiv
!>          (sqrt(2/3 e':e'), always of the deviatoric part: v2's defrate_eq receives the
!>          deviatoric d from velgrad_sympart, and strain_equiv receives deviatoric 5-vectors).
!>
!> Representations (all real(rk)):
!>   Mandel 6-vector of a symmetric A:  {a} = {A11, A22, A33, sqrt2 A23, sqrt2 A13, sqrt2 A12}.
!>     The Mandel basis is orthonormal, so A:B = dot_product({a}, {b}) and |A| = |{a}|.
!>     Component p corresponds to the tensor indices (mandel_index_i(p), mandel_index_j(p)).
!>   Tensor-component 6-vector (I/O only, conventions §1.2): {A11, A22, A33, A23, A13, A12}.
!>   Skew 3-vector of a skew W: {w} = {W23, W13, W12}, i.e.
!>     W = [[0, w12, w13], [-w12, 0, w23], [-w13, -w23, 0]]. No axial-vector sign is implied.
!>   Deviatoric 5-vector of A (addendum_5vectors.md (F1), KD25 Eq. 8):
!>     {a'} = {(a11 - a22)/sqrt2, sqrt(3/2) a'33, sqrt2 a23, sqrt2 a13, sqrt2 a12}, the
!>     components of A' in an orthonormal basis of deviatoric tensors: A':B' = {a'}.{b'}.
!>   6x6 Mandel matrices act on Mandel 6-vectors; the inner product (dot_product) of two Mandel
!>   vectors and the product of a 6x6 Mandel matrix with a Mandel vector (matmul) need no
!>   routine of their own.
!>
!> All routines are pure functions of one tensor (one material point), for use inside loops
!> over points (conventions §6). None depends on a frame: the caller's variables carry the
!> frame suffixes.
!>
!> Public parameters:
!>   mandel_index_i  integer(6)  first tensor index of Mandel component p:  1 2 3 2 1 1
!>   mandel_index_j  integer(6)  second tensor index of Mandel component p: 1 2 3 3 3 2
!>   dev5_projection real(rk)(5,6)  [P]: {a'} = [P]{a} (addendum (F3)); [P][P]^T = I5 and
!>                                  [P]^T[P] is the deviatoric projector
!>
!> Public routines:
!>   mandel_from_matrix      3x3 matrix (its symmetric part) -> Mandel 6-vector
!>   matrix_from_mandel      Mandel 6-vector -> symmetric 3x3 matrix
!>   mandel_from_components  tensor-component 6-vector -> Mandel 6-vector
!>   components_from_mandel  Mandel 6-vector -> tensor-component 6-vector
!>   mandel_trace            trace of the tensor, A11 + A22 + A33
!>   mandel_dev              deviatoric part, A - (tr A / 3) I
!>   mandel_norm             Frobenius norm, sqrt(A:A)
!>   equiv_stress            von Mises equivalent stress, sqrt(3/2 s':s')
!>   equiv_strain            equivalent strain (or strain rate), sqrt(2/3 e':e')
!>   skew_from_matrix        3x3 matrix (its skew part) -> skew 3-vector
!>   matrix_from_skew        skew 3-vector -> skew 3x3 matrix
!>   mandel_rotation         [Q6(R)]: {R A R^T} = [Q6(R)] {A}, Eq. (R4)
!>   spin_commutator         [W^]: {A W - W A} = [W^] {A}, Eq. (R8)
!>   dev5_from_mandel        Mandel 6-vector -> deviatoric 5-vector, [P]{a}
!>   mandel_from_dev5        deviatoric 5-vector -> Mandel 6-vector of A', [P]^T{a'}
!>   mat5_from_mat6          5x5 form [P][M][P]^T of a 6x6 Mandel matrix
module tensor_mod

  use kinds_mod, only: rk
  use constants_mod, only: sqrt2, inv_sqrt2, inv_sqrt6, sqrt_3_2, sqrt_2_3, one_third, &
    & identity_mandel

  implicit none

  private

  public :: mandel_index_i, mandel_index_j
  public :: mandel_from_matrix, matrix_from_mandel
  public :: mandel_from_components, components_from_mandel
  public :: mandel_trace, mandel_dev, mandel_norm
  public :: equiv_stress, equiv_strain
  public :: skew_from_matrix, matrix_from_skew
  public :: mandel_rotation, spin_commutator
  public :: dev5_projection
  public :: dev5_from_mandel, mandel_from_dev5, mat5_from_mat6

  ! Tensor indices (i, j) of the Mandel components, order 11 22 33 23 13 12.
  integer, parameter :: mandel_index_i(6) = [1, 2, 3, 2, 1, 1]
  integer, parameter :: mandel_index_j(6) = [1, 2, 3, 3, 3, 2]

  ! Scale factor of each Mandel component relative to the tensor component: 1 on the normal
  ! components, sqrt(2) on the shear components.
  real(rk), parameter :: mandel_scale(6) = [1.0_rk, 1.0_rk, 1.0_rk, sqrt2, sqrt2, sqrt2]

  ! Projection [P] from Mandel 6-vectors to deviatoric 5-vectors (addendum_5vectors.md (F3);
  ! Mandel form of KD25 Eq. 28). Its rows are the Mandel vectors of the deviatoric basis
  ! B1..B5. Stored column by column (reshape fills columns first).
  real(rk), parameter :: dev5_projection(5, 6) = reshape([ &
    &  inv_sqrt2, -inv_sqrt6, 0.0_rk, 0.0_rk, 0.0_rk, &
    & -inv_sqrt2, -inv_sqrt6, 0.0_rk, 0.0_rk, 0.0_rk, &
    &  0.0_rk, 2.0_rk*inv_sqrt6, 0.0_rk, 0.0_rk, 0.0_rk, &
    &  0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, &
    &  0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, &
    &  0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk], [5, 6])

contains

  !> Mandel 6-vector of the symmetric part of a 3x3 matrix.
  !>
  !> Theory : conventions.md §1.1. The shear components use the average of the two
  !>          off-diagonal entries, so a non-symmetric input gives its symmetric part.
  !>
  !> Arguments:
  !>   mat  (in)  real(rk)(3,3)  matrix (normally symmetric)
  !>   vec  (out, result) real(rk)(6)  Mandel 6-vector of sym(mat)
  pure function mandel_from_matrix(mat) result(vec)

    real(rk), intent(in) :: mat(3, 3)
    real(rk) :: vec(6)

    ! Normal components unchanged.
    vec(1) = mat(1, 1)
    vec(2) = mat(2, 2)
    vec(3) = mat(3, 3)

    ! Shear components: sqrt(2) times the symmetric part, sqrt(2) (A_ij + A_ji)/2.
    vec(4) = inv_sqrt2*(mat(2, 3) + mat(3, 2))
    vec(5) = inv_sqrt2*(mat(1, 3) + mat(3, 1))
    vec(6) = inv_sqrt2*(mat(1, 2) + mat(2, 1))

  end function mandel_from_matrix

  !> Symmetric 3x3 matrix of a Mandel 6-vector.
  !>
  !> Theory : conventions.md §1.1.
  !>
  !> Arguments:
  !>   vec  (in)  real(rk)(6)    Mandel 6-vector
  !>   mat  (out, result) real(rk)(3,3)  symmetric matrix
  pure function matrix_from_mandel(vec) result(mat)

    real(rk), intent(in) :: vec(6)
    real(rk) :: mat(3, 3)

    ! Diagonal.
    mat(1, 1) = vec(1)
    mat(2, 2) = vec(2)
    mat(3, 3) = vec(3)

    ! Off-diagonal: tensor component = Mandel component / sqrt(2), both triangles.
    mat(2, 3) = inv_sqrt2*vec(4)
    mat(1, 3) = inv_sqrt2*vec(5)
    mat(1, 2) = inv_sqrt2*vec(6)
    mat(3, 2) = mat(2, 3)
    mat(3, 1) = mat(1, 3)
    mat(2, 1) = mat(1, 2)

  end function matrix_from_mandel

  !> Mandel 6-vector of a tensor-component 6-vector (input/output form, conventions §1.2).
  !>
  !> Theory : conventions.md §1.1–1.2.
  !>
  !> Arguments:
  !>   comp  (in)  real(rk)(6)  tensor components A11 A22 A33 A23 A13 A12
  !>   vec   (out, result) real(rk)(6)  Mandel 6-vector
  pure function mandel_from_components(comp) result(vec)

    real(rk), intent(in) :: comp(6)
    real(rk) :: vec(6)

    ! Scale the shear components by sqrt(2).
    vec(1:3) = comp(1:3)
    vec(4:6) = sqrt2*comp(4:6)

  end function mandel_from_components

  !> Tensor-component 6-vector (input/output form, conventions §1.2) of a Mandel 6-vector.
  !>
  !> Theory : conventions.md §1.1–1.2.
  !>
  !> Arguments:
  !>   vec   (in)  real(rk)(6)  Mandel 6-vector
  !>   comp  (out, result) real(rk)(6)  tensor components A11 A22 A33 A23 A13 A12
  pure function components_from_mandel(vec) result(comp)

    real(rk), intent(in) :: vec(6)
    real(rk) :: comp(6)

    ! Divide the shear components by sqrt(2).
    comp(1:3) = vec(1:3)
    comp(4:6) = inv_sqrt2*vec(4:6)

  end function components_from_mandel

  !> Trace of the tensor represented by a Mandel 6-vector.
  !>
  !> Theory : conventions.md §1.1 (tr A = A11 + A22 + A33 = I:A).
  !>
  !> Arguments:
  !>   vec    (in)  real(rk)(6)  Mandel 6-vector
  !>   trace  (out, result) real(rk)  trace of the tensor
  pure function mandel_trace(vec) result(trace)

    real(rk), intent(in) :: vec(6)
    real(rk) :: trace

    ! Sum of the normal components.
    trace = vec(1) + vec(2) + vec(3)

  end function mandel_trace

  !> Deviatoric part of the tensor represented by a Mandel 6-vector, A' = A - (tr A / 3) I.
  !>
  !> Theory : conventions.md §1.1; the Mandel vector of I is {1,1,1,0,0,0}.
  !>
  !> Arguments:
  !>   vec      (in)  real(rk)(6)  Mandel 6-vector of A
  !>   vec_dev  (out, result) real(rk)(6)  Mandel 6-vector of A'
  pure function mandel_dev(vec) result(vec_dev)

    real(rk), intent(in) :: vec(6)
    real(rk) :: vec_dev(6)

    ! Subtract the mean normal component from the normal components.
    vec_dev = vec - (one_third*mandel_trace(vec))*identity_mandel

  end function mandel_dev

  !> Frobenius norm of the tensor represented by a Mandel 6-vector, |A| = sqrt(A:A).
  !>
  !> Theory : conventions.md §1.1 (the Mandel basis is orthonormal).
  !>
  !> Arguments:
  !>   vec   (in)  real(rk)(6)  Mandel 6-vector
  !>   norm  (out, result) real(rk)  Frobenius norm of the tensor
  pure function mandel_norm(vec) result(norm)

    real(rk), intent(in) :: vec(6)
    real(rk) :: norm

    ! Euclidean norm of the Mandel vector.
    norm = sqrt(dot_product(vec, vec))

  end function mandel_norm

  !> Von Mises equivalent stress, sqrt(3/2 s':s'), s' the deviatoric stress.
  !>
  !> Source : v2 matrix_operations_mod.f90::stress_equiv_3x3 (same definition).
  !>
  !> Arguments:
  !>   stress  (in)  real(rk)(6)  stress, Mandel 6-vector [MPa]
  !>   eq      (out, result) real(rk)  equivalent stress [MPa]
  pure function equiv_stress(stress) result(eq)

    real(rk), intent(in) :: stress(6)
    real(rk) :: eq

    ! sqrt(3/2) times the norm of the deviatoric part.
    eq = sqrt_3_2*mandel_norm(mandel_dev(stress))

  end function equiv_stress

  !> Equivalent strain or strain rate, sqrt(2/3 e':e'), e' the deviatoric part.
  !>
  !> Source : v2 matrix_operations_mod.f90::strain_equiv_3x3 (strain, elastic and plastic
  !>          strain), ::defrate_eq (deformation rate; v2 passes the deviatoric part) and
  !>          ::strain_equiv (plastic deformation rate as a deviatoric 5-vector): same definition.
  !>
  !> Arguments:
  !>   strain  (in)  real(rk)(6)  strain [-] or strain rate [1/s], Mandel 6-vector
  !>   eq      (out, result) real(rk)  equivalent strain [-] or strain rate [1/s]
  pure function equiv_strain(strain) result(eq)

    real(rk), intent(in) :: strain(6)
    real(rk) :: eq

    ! sqrt(2/3) times the norm of the deviatoric part.
    eq = sqrt_2_3*mandel_norm(mandel_dev(strain))

  end function equiv_strain

  !> Skew 3-vector {W23, W13, W12} of the skew part W = (A - A^T)/2 of a 3x3 matrix.
  !>
  !> Theory : conventions.md §1.4.
  !>
  !> Arguments:
  !>   mat  (in)  real(rk)(3,3)  matrix
  !>   skw  (out, result) real(rk)(3)  skew 3-vector of skew(mat)
  pure function skew_from_matrix(mat) result(skw)

    real(rk), intent(in) :: mat(3, 3)
    real(rk) :: skw(3)

    ! Upper-triangle components of the skew part, order 23 13 12.
    skw(1) = 0.5_rk*(mat(2, 3) - mat(3, 2))
    skw(2) = 0.5_rk*(mat(1, 3) - mat(3, 1))
    skw(3) = 0.5_rk*(mat(1, 2) - mat(2, 1))

  end function skew_from_matrix

  !> Skew 3x3 matrix of a skew 3-vector {w23, w13, w12}.
  !>
  !> Theory : conventions.md §1.4, W = [[0, w12, w13], [-w12, 0, w23], [-w13, -w23, 0]].
  !>
  !> Arguments:
  !>   skw  (in)  real(rk)(3)    skew 3-vector
  !>   mat  (out, result) real(rk)(3,3)  skew matrix
  pure function matrix_from_skew(skw) result(mat)

    real(rk), intent(in) :: skw(3)
    real(rk) :: mat(3, 3)

    ! Zero diagonal.
    mat(1, 1) = 0.0_rk
    mat(2, 2) = 0.0_rk
    mat(3, 3) = 0.0_rk

    ! Upper triangle from the vector, lower triangle by antisymmetry.
    mat(2, 3) = skw(1)
    mat(1, 3) = skw(2)
    mat(1, 2) = skw(3)
    mat(3, 2) = -skw(1)
    mat(3, 1) = -skw(2)
    mat(2, 1) = -skw(3)

  end function matrix_from_skew

  !> Mandel rotation operator [Q6(R)] of a 3x3 rotation R: {R A R^T} = [Q6(R)] {A} for every
  !> symmetric A. [Q6] is orthogonal when R is, and [Q6(R)]^T = [Q6(R^T)].
  !>
  !> Theory : recoupled_formulation.md Eq. (R4): column q is the Mandel vector of R E_q R^T,
  !>          E_q the q-th Mandel basis tensor (E_q = e_i (x) e_i for q <= 3,
  !>          (e_i (x) e_j + e_j (x) e_i)/sqrt(2) for q >= 4, (i, j) = (mandel_index_i(q),
  !>          mandel_index_j(q))). Writing (k, l) for the indices of row p, this gives in closed
  !>          form
  !>            Q6(p, q) = (s_p s_q / 2) (R_ki R_lj + R_kj R_li),
  !>          s_p = mandel_scale(p) (1 normal, sqrt(2) shear) (Mehrabadi & Cowin 1990). The
  !>          closed form is checked against the column definition by the unit tests.
  !>
  !> Arguments:
  !>   rot  (in)  real(rk)(3,3)  rotation matrix R (orthogonal, det +1; not checked)
  !>   q6   (out, result) real(rk)(6,6)  Mandel rotation operator [Q6(R)]
  !> Locals:
  !>   p, q   integer  row and column of [Q6]
  !>   i, j   integer  tensor indices of Mandel component q
  !>   k, l   integer  tensor indices of Mandel component p
  pure function mandel_rotation(rot) result(q6)

    real(rk), intent(in) :: rot(3, 3)
    real(rk) :: q6(6, 6)

    integer :: p, q
    integer :: i, j, k, l

    ! Fill [Q6] entry by entry from the closed form.
    do q = 1, 6
      i = mandel_index_i(q)
      j = mandel_index_j(q)
      do p = 1, 6
        k = mandel_index_i(p)
        l = mandel_index_j(p)
        q6(p, q) = 0.5_rk*mandel_scale(p)*mandel_scale(q) &
          & *(rot(k, i)*rot(l, j) + rot(k, j)*rot(l, i))
      end do
    end do

  end function mandel_rotation

  !> Spin-commutator matrix [W^] of a skew tensor W: {A W - W A} = [W^] {A} for every
  !> symmetric A (Mandel vectors). [W^] is antisymmetric, and its output is traceless.
  !>
  !> Theory : recoupled_formulation.md Eq. (R8), used there with A = e^e and W = w^p:
  !>          {e^e w^p - w^p e^e} = [W^]{e^e}. Rows/columns in Mandel order 11 22 33 23 13 12.
  !>
  !> Arguments:
  !>   skw   (in)  real(rk)(3)    skew 3-vector of W, {w23, w13, w12} (conventions §1.4)
  !>   spin  (out, result) real(rk)(6,6)  spin-commutator matrix [W^]
  !> Locals:
  !>   w23, w13, w12   real(rk)  components of W
  !>   s23, s13, s12   real(rk)  sqrt(2) times the components of W
  pure function spin_commutator(skw) result(spin)

    real(rk), intent(in) :: skw(3)
    real(rk) :: spin(6, 6)

    real(rk) :: w23, w13, w12
    real(rk) :: s23, s13, s12

    ! Components of W and their sqrt(2) multiples.
    w23 = skw(1)
    w13 = skw(2)
    w12 = skw(3)
    s23 = sqrt2*w23
    s13 = sqrt2*w13
    s12 = sqrt2*w12

    ! Eq. (R8), written row by row (spin(row, :)).
    spin(1, :) = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, -s13, -s12]
    spin(2, :) = [0.0_rk, 0.0_rk, 0.0_rk, -s23, 0.0_rk, s12]
    spin(3, :) = [0.0_rk, 0.0_rk, 0.0_rk, s23, s13, 0.0_rk]
    spin(4, :) = [0.0_rk, s23, -s23, 0.0_rk, w12, w13]
    spin(5, :) = [s13, 0.0_rk, -s13, -w12, 0.0_rk, w23]
    spin(6, :) = [s12, -s12, 0.0_rk, -w13, -w23, 0.0_rk]

  end function spin_commutator

  !> Deviatoric 5-vector of the tensor represented by a Mandel 6-vector, {a'} = [P]{a}. The
  !> spherical part of the input is discarded (no prior deviatoric projection is needed).
  !>
  !> Theory : addendum_5vectors.md (F1), (F3); KD25 Eqs. 8 and 28.
  !>
  !> Arguments:
  !>   vec   (in)  real(rk)(6)  Mandel 6-vector of A
  !>   dev5  (out, result) real(rk)(5)  deviatoric 5-vector of A'
  pure function dev5_from_mandel(vec) result(dev5)

    real(rk), intent(in) :: vec(6)
    real(rk) :: dev5(5)

    ! (a11 - a22)/sqrt(2).
    dev5(1) = inv_sqrt2*(vec(1) - vec(2))

    ! sqrt(3/2) a'33 = (2 a33 - a11 - a22)/sqrt(6).
    dev5(2) = inv_sqrt6*(2.0_rk*vec(3) - vec(1) - vec(2))

    ! Shear components: the Mandel components (sqrt(2) a23, sqrt(2) a13, sqrt(2) a12).
    dev5(3:5) = vec(4:6)

  end function dev5_from_mandel

  !> Mandel 6-vector of the deviatoric tensor A' represented by a 5-vector, [P]^T{a'}.
  !>
  !> Theory : addendum_5vectors.md (F4); KD25 Eq. 29 (with the Mandel scaling of the shears).
  !>
  !> Arguments:
  !>   dev5  (in)  real(rk)(5)  deviatoric 5-vector
  !>   vec   (out, result) real(rk)(6)  Mandel 6-vector of A' (traceless)
  pure function mandel_from_dev5(dev5) result(vec)

    real(rk), intent(in) :: dev5(5)
    real(rk) :: vec(6)

    ! Normal components: columns 1-2 of [P]^T.
    vec(1) = inv_sqrt2*dev5(1) - inv_sqrt6*dev5(2)
    vec(2) = -inv_sqrt2*dev5(1) - inv_sqrt6*dev5(2)
    vec(3) = 2.0_rk*inv_sqrt6*dev5(2)

    ! Shear components unchanged.
    vec(4:6) = dev5(3:5)

  end function mandel_from_dev5

  !> 5x5 form of a 6x6 Mandel matrix, [M'] = [P][M][P]^T: its deviatoric-deviatoric block in
  !> the 5-vector basis. Exact for operators that preserve the deviatoric and spherical
  !> subspaces ([Q6], [W^], [M*]); for an anisotropic stiffness it is v2's decoupled form.
  !>
  !> Theory : addendum_5vectors.md §4, (F7) (KD25 Eq. 30 generalized).
  !>
  !> Arguments:
  !>   mat6  (in)  real(rk)(6,6)  Mandel matrix [M]
  !>   mat5  (out, result) real(rk)(5,5)  5x5 form [M']
  pure function mat5_from_mat6(mat6) result(mat5)

    real(rk), intent(in) :: mat6(6, 6)
    real(rk) :: mat5(5, 5)

    ! [P][M][P]^T.
    mat5 = matmul(dev5_projection, matmul(mat6, transpose(dev5_projection)))

  end function mat5_from_mat6

end module tensor_mod
