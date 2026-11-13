! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Vertices of the single crystal yield surface, which start the increment-1 viscoplastic
!> estimate (crystal_plasticity_model_mod::cp_vertex_estimate): computed from a phase's slip
!> systems and strengths (initialization isotropic_viscoplastic, ROADMAP step 6.7, D57, D78), as
!> deviatoric Mandel 6-vectors (order 11 22 33 23 13 12), crystal frame. (FEPX v2's tabulated
!> vertices, kept for the comparison with v2, were removed in step 9.4, D56.)
!>
!> Theory : dev/theory/addendum_yield_vertices.md (approved 2026-10-11), after Kasemer, M. and
!>          Dawson, P.R. (2025), arXiv:2502.19531, refs/kasemer_dawson_2025.pdf (KD25), §3,
!>          Eqs. 6-13 and Algorithm 1. Schmid's law |p^a . s| = tau_0^a is a pair of hyperplanes
!>          in the deviatoric 5-vector space (KD25 Eqs. 6-9; addendum_5vectors.md (F1)); every
!>          set of 5 linearly independent systems meets in a candidate vertex for each sign
!>          pattern of the strengths (KD25 Eqs. 10-12); a candidate inside or on every
!>          hyperplane is a vertex of the yield surface (KD25 Eq. 13, the inner convex hull).
!>          The strengths enter KD25 Eq. 12 directly (user, D57): the vertices are stresses.
!>          Strengths are the same in both slip directions in FEPX, so the vertices come in
!>          pairs +-s; one of each pair is kept (as v2's tables), the caller applies the sign.
!>
!> Public parameters (tolerances of addendum §3, each 5-6 orders of magnitude from both sides of
!> the gaps measured on FEPX's slip systems):
!>   yield_tol_det   1e-10  |det P| relative to Hadamard's bound prod |p^a| (T1): measured
!>                          dependent sets <= 1.4e-16, independent >= 3.1e-3 (bct)
!>   yield_tol_hull  1e-10  |p^b . s| - tau_0^b relative to tau_0^b, the hull test and the residual
!>                          of the 5 defining equations (T2): measured vertices <= 5.3e-15,
!>                          rejected candidates >= 5.4e-3 (bct)
!>   yield_tol_dup   1e-8   distance of two vertices relative to the smallest strength (T3):
!>                          measured copies <= 1.4e-14 apart, distinct vertices >= 3.3e-2 (bct)
!>
!> Public routines:
!>   yield_vertices_compute  the vertices (one of each +- pair) of a set of slip systems [MPa]
module yield_vertices_mod

  use kinds_mod, only: rk
  use tensor_mod, only: dev5_from_mandel, mandel_from_dev5
  use linalg_mod, only: lu_factor, lu_solve

  implicit none

  private

  public :: yield_vertices_compute
  public :: yield_tol_det, yield_tol_hull, yield_tol_dup

  real(rk), parameter :: yield_tol_det = 1.0e-10_rk
  real(rk), parameter :: yield_tol_hull = 1.0e-10_rk
  real(rk), parameter :: yield_tol_dup = 1.0e-8_rk

contains

  !> The vertices of the single crystal yield surface of N slip systems with strengths tau_0^a
  !> (the same in both slip directions), one of each +- pair, as deviatoric Mandel 6-vectors.
  !>
  !> Theory : addendum_yield_vertices.md §2 (steps V-1 to V-7), KD25 Algorithm 1:
  !>          V-1 p^a = 5-vector of the Schmid tensor; V-3 for each combination a<b<c<d<e
  !>          (lexicographic), [P] = rows p^a..p^e, independent when
  !>          |det P| > yield_tol_det prod |p|; V-4 for the 16 sign patterns with + in the first
  !>          row, [P]{s} = {+-tau} (LU with partial pivoting, one factorization), accepted when
  !>          the five equations hold to yield_tol_hull tau; V-5 a vertex when
  !>          |p^b . s| - tau_0^b <= yield_tol_hull tau_0^b for every b (stops at the first
  !>          violation); V-6 kept unless it or its negative lies within yield_tol_dup tau_0_min of
  !>          a vertex already kept (first found is kept). V-7 (fewer than 5 independent
  !>          systems: no vertex) is reported through num_independent = 0; the caller decides.
  !> Source : KD25 §3 (method); Golub & Van Loan (2013) Alg. 3.4.1 through linalg_mod (solve).
  !>
  !> Arguments:
  !>   schmid_sym_c     (in)  real(rk)(6, N)    Schmid tensors p of the systems, Mandel, crystal
  !>                                            frame (crystal_mod) [-]
  !>   strength         (in)  real(rk)(N)       strengths tau_0^a > 0 [MPa]
  !>   vertices_c       (out) real(rk)(6, :)    the vertices, one of each +- pair, deviatoric
  !>                                            Mandel, crystal frame, in the order found [MPa]
  !>   num_independent  (out) integer           number of linearly independent combinations
  !>                                            M_L (0: the systems span fewer than 5
  !>                                            dimensions, no vertex)
  !> Locals:
  !>   n             integer               number of systems N
  !>   p5            real(rk)(5, N)        Schmid 5-vectors p^a
  !>   p_norm        real(rk)(N)           |p^a|
  !>   tau_0_min     real(rk)              smallest strength [MPa]
  !>   combo         integer(5)            current combination of systems
  !>   mat_p, lu     real(rk)(5, 5)        [P] and its LU factors
  !>   perm          integer(5)            row permutation of the factorization
  !>   det           real(rk)              det P
  !>   rhs, s        real(rk)(5)           signed strengths, candidate vertex [MPa]
  !>   i_pattern     integer               sign pattern 0..15 (bit k-2: sign of row k, k >= 2)
  !>   k, b, status  integer               row, system, factorization status
  !>   is_vertex     logical               the candidate passed the residual and hull tests
  !>   kept          real(rk)(5, :)        the vertices kept so far (5-vectors) [MPa]
  !>   num_kept      integer               their number
  !>   more          logical               another combination exists
  subroutine yield_vertices_compute(schmid_sym_c, strength, vertices_c, num_independent)

    real(rk), intent(in) :: schmid_sym_c(:, :)
    real(rk), intent(in) :: strength(:)
    real(rk), allocatable, intent(out) :: vertices_c(:, :)
    integer, intent(out) :: num_independent

    integer :: n, combo(5), perm(5), i_pattern, k, b, status, num_kept
    real(rk) :: p5(5, size(strength)), p_norm(size(strength)), tau_0_min
    real(rk) :: mat_p(5, 5), lu(5, 5), det, rhs(5), s(5)
    real(rk), allocatable :: kept(:, :)
    logical :: is_vertex, more

    n = size(strength)
    num_independent = 0
    num_kept = 0
    allocate (kept(5, 64))

    ! V-1: the Schmid 5-vectors and their norms (1/sqrt(2) for unit s and m).
    do b = 1, n
      p5(:, b) = dev5_from_mandel(schmid_sym_c(:, b))
      p_norm(b) = norm2(p5(:, b))
    end do
    tau_0_min = minval(strength)

    ! V-3: every combination of 5 systems, in lexicographic order (none if N < 5).
    combo = [1, 2, 3, 4, 5]
    more = n >= 5
    do while (more)
      ! [P] with rows p^a..p^e, its factors and determinant; skip a dependent set (V2).
      do k = 1, 5
        mat_p(k, :) = p5(:, combo(k))
      end do
      call lu_factor(mat_p, lu, perm, det, status)
      if (status == 0 .and. abs(det) > yield_tol_det*product(p_norm(combo))) then
        num_independent = num_independent + 1

        ! V-4: the 16 sign patterns with +tau_0 in the first row.
        do i_pattern = 0, 15
          rhs(1) = strength(combo(1))
          do k = 2, 5
            if (btest(i_pattern, k - 2)) then
              rhs(k) = -strength(combo(k))
            else
              rhs(k) = strength(combo(k))
            end if
          end do
          call lu_solve(lu, perm, rhs, s)

          ! The five defining equations to yield_tol_hull (rejects an ill-conditioned solve).
          is_vertex = all(abs(matmul(mat_p, s) - rhs) <= yield_tol_hull*strength(combo))

          ! V-5: inside or on every hyperplane pair, stopping at the first violation.
          if (is_vertex) then
            do b = 1, n
              if (abs(dot_product(p5(:, b), s)) - strength(b) > yield_tol_hull*strength(b)) then
                is_vertex = .false.
                exit
              end if
            end do
          end if

          ! V-6: keep it unless it or its negative is already kept.
          if (is_vertex) call keep_vertex(s)
        end do
      end if

      call next_combination(combo, n, more)
    end do

    ! Output: deviatoric Mandel 6-vectors.
    allocate (vertices_c(6, num_kept))
    do k = 1, num_kept
      vertices_c(:, k) = mandel_from_dev5(kept(:, k))
    end do

  contains

    !> Append a vertex to the kept list unless it or its negative lies within
    !> yield_tol_dup tau_0_min of a kept vertex (the list grows by doubling).
    !>
    !> Arguments:
    !>   vertex  (in)  real(rk)(5)  candidate vertex, 5-vector [MPa]
    !> Locals:
    !>   j       integer               kept vertex
    !>   grown   real(rk)(5, :)        the list with doubled capacity
    subroutine keep_vertex(vertex)

      real(rk), intent(in) :: vertex(5)

      integer :: j
      real(rk), allocatable :: grown(:, :)

      ! A copy of a kept vertex, or of its negative: nothing to do.
      do j = 1, num_kept
        if (norm2(kept(:, j) - vertex) <= yield_tol_dup*tau_0_min .or. &
          & norm2(kept(:, j) + vertex) <= yield_tol_dup*tau_0_min) return
      end do

      ! Append, doubling the capacity when full.
      if (num_kept == size(kept, 2)) then
        allocate (grown(5, 2*size(kept, 2)))
        grown(:, 1:num_kept) = kept(:, 1:num_kept)
        call move_alloc(grown, kept)
      end if
      num_kept = num_kept + 1
      kept(:, num_kept) = vertex

    end subroutine keep_vertex

  end subroutine yield_vertices_compute

  !> The next combination of 5 out of 1..n in lexicographic order.
  !>
  !> Source : the standard successor of a k-combination (rightmost entry that can be increased
  !>          is increased, the following entries set consecutively); no numerical method.
  !>
  !> Arguments:
  !>   combo  (inout) integer(5)  current combination, increasing entries in 1..n
  !>   n      (in)    integer     number of systems
  !>   more   (out)   logical     .false. when combo was the last combination (combo unchanged)
  !> Locals:
  !>   k, j   integer             entry increased, following entries
  pure subroutine next_combination(combo, n, more)

    integer, intent(inout) :: combo(5)
    integer, intent(in) :: n
    logical, intent(out) :: more

    integer :: k, j

    ! The rightmost entry k that is below its largest value n - 5 + k.
    more = .false.
    do k = 5, 1, -1
      if (combo(k) < n - 5 + k) then
        ! Increase it and set the following entries consecutively.
        combo(k) = combo(k) + 1
        do j = k + 1, 5
          combo(j) = combo(j - 1) + 1
        end do
        more = .true.
        return
      end if
    end do

  end subroutine next_combination

end module yield_vertices_mod
