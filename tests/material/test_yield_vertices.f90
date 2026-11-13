! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of the computed single crystal yield-surface vertices
!> (yield_vertices_mod::yield_vertices_compute; initialization isotropic_viscoplastic, ROADMAP step 6.7, D57).
!>
!> Source : FEPX (ROADMAP step 6.7; dev/theory/addendum_yield_vertices.md §7).
!>
!> Checks:
!>   - against KD25's own implementation (fixture material/fixtures/kd25_vertices.txt, README):
!>     fcc unit strengths (KD25 §5.1: 56 vertices = 28 pairs), bcc {110}<111> with system 1
!>     5 % weaker (KD25 §6.2: 116 = 58 pairs), hcp c/a 1.587 with 1 / 1.25 / 1.7, hcp c/a 1.59
!>     pyramidal <c+a> only through the active families (KD25 §5.2: 92 = 46 pairs), bct c/a 1:
!>     the same number of vertices and the same set (+- pairs, any order) to 1e-3 of the
!>     strength (KD25's code rounds to 4 decimals);
!>   - fcc at unit strength: 28 pairs from 384 independent combinations (the size of v2's fcc
!>     table, which matched the computed set to its rounding, 2e-5, until the table was removed
!>     in step 9.4);
!>   - properties for fcc, bcc 24 systems, hcp c/a 1.587 and bct c/a 0.5456 with random
!>     strengths 100-200 MPa per system (fixed seed): every vertex deviatoric (Mandel trace
!>     0) and inside every hyperplane pair, |p^b : s| <= tau^b (1 + 1e-10), with p^b : s
!>     evaluated by explicit 3x3 algebra (the slip direction and normal, d . S . n) from the
!>     Mandel 6-vector, so the stored vertices are Mandel 6-vectors of the stresses; the
!>     systems at their strength (to 1e-9) span the 5-dimensional deviatoric space (det of
!>     their 5-vector Gram matrix relative to its scale > 1e-12); no two vertices, or one and
!>     the negative of another, closer than yield_tol_dup tau_min; scaling every strength by 3
!>     scales the vertices by 3 (relative 1e-12, the same order); reversing the order of the
!>     systems gives the same set (1e-10 relative);
!>   - fewer than 5 independent systems (hcp basal + prismatic): no independent combination,
!>     no vertex.
!>
!> Locals:
!>   fixture  character  path of kd25_vertices.txt
!>   crys     type(crystal_t)  slip geometry
program test_yield_vertices

  use kinds_mod, only: rk
  use tensor_mod, only: mandel_trace, matrix_from_mandel, dev5_from_mandel
  use linalg_mod, only: lu_factor
  use crystal_mod, only: crystal_t, crystal_init, crystal_fcc, crystal_bcc, crystal_hcp, &
    & crystal_bct
  use yield_vertices_mod, only: yield_vertices_compute, yield_tol_dup
  use tensor_test_utils_mod, only: seed_random
  use assert_mod, only: assert_true, assert_equal, test_summary

  implicit none

  character(len=1024) :: fixture
  type(crystal_t) :: crys

  call get_command_argument(1, fixture)
  call seed_random(20261011)

  call check_fixture()
  call check_fcc_counts()

  call crystal_init(crys, crystal_fcc)
  call check_properties('fcc')
  call crystal_init(crys, crystal_bcc, families=[1, 2])
  call check_properties('bcc24')
  call crystal_init(crys, crystal_hcp, c_over_a=1.587_rk)
  call check_properties('hcp')
  call crystal_init(crys, crystal_bct, c_over_a=0.5456_rk)
  call check_properties('bct')

  call check_rank_deficient()

  call test_summary()

contains

  !> The five KD25 cases of the fixture (header).
  !>
  !> Locals:
  !>   unit, ios      integer          file unit, status
  !>   structure      character(3)     crystal structure of a case
  !>   c_over_a       real(rk)         c/a of a case
  !>   num_slip, num_ref  integer      systems, KD25's number of vertices
  !>   strength       real(rk)(:)      strengths of the case
  !>   ref5           real(rk)(:, :)   KD25's vertices, 5-vectors (both signs)
  !>   vertices       real(rk)(:, :)   computed vertices, Mandel (one of each pair) [MPa]
  !>   num_indep      integer          independent combinations
  !>   i_case, i      integer          case, vertex
  !>   label          character        case label
  subroutine check_fixture()

    integer :: unit, ios, num_slip, num_ref, num_indep, i_case, i
    character(len=3) :: structure
    character(len=32) :: label
    real(rk) :: c_over_a
    real(rk), allocatable :: strength(:), ref5(:, :), vertices(:, :)

    open (newunit=unit, file=trim(fixture), status='old', action='read', iostat=ios)
    call assert_equal(ios, 0, 'fixture opened')
    if (ios /= 0) return

    do i_case = 1, 5
      read (unit, *) structure, c_over_a, num_slip, num_ref
      allocate (strength(num_slip), ref5(5, num_ref))
      read (unit, *) strength
      do i = 1, num_ref
        read (unit, *) ref5(:, i)
      end do
      write (label, '(a,a,i0)') trim(structure), ' case ', i_case

      ! The slip systems of the case (pyramidal <c+a> only for case 4).
      select case (i_case)
      case (1)
        call crystal_init(crys, crystal_fcc)
      case (2)
        call crystal_init(crys, crystal_bcc)
      case (3)
        call crystal_init(crys, crystal_hcp, c_over_a=c_over_a)
      case (4)
        call crystal_init(crys, crystal_hcp, c_over_a=c_over_a, families=[3])
        strength = strength(7:18)
      case default
        call crystal_init(crys, crystal_bct, c_over_a=c_over_a)
      end select

      call yield_vertices_compute(crys%schmid_sym_c, strength, vertices, num_indep)
      call assert_equal(2*size(vertices, 2), num_ref, trim(label) // ': number of vertices')
      call assert_true(same_set(vertices, mandel_from_ref(ref5), 1.0e-3_rk*minval(strength)), &
        & trim(label) // ': the vertices of KD25''s implementation')
      deallocate (strength, ref5)
    end do
    close (unit)

  end subroutine check_fixture

  !> The fcc vertex and combination counts at unit strength (header).
  !>
  !> Locals:
  !>   vertices   real(rk)(:, :)  computed vertices
  !>   num_indep  integer         independent combinations
  subroutine check_fcc_counts()

    real(rk), allocatable :: vertices(:, :)
    integer :: num_indep

    call crystal_init(crys, crystal_fcc)
    call yield_vertices_compute(crys%schmid_sym_c, spread(1.0_rk, 1, 12), vertices, num_indep)
    call assert_equal(size(vertices, 2), 28, 'fcc: 28 pairs')
    call assert_equal(num_indep, 384, 'fcc: 384 independent combinations')

  end subroutine check_fcc_counts

  !> Properties of the vertices of the current crystal with random strengths (header).
  !>
  !> Arguments:
  !>   name  (in)  character  label
  !> Locals:
  !>   n                integer         systems
  !>   strength         real(rk)(n)     random strengths 100-200 [MPa]
  !>   vertices, scaled, reversed  real(rk)(:, :)  vertices of the strengths, of 3 x the
  !>                                    strengths, of the reversed systems [MPa]
  !>   num_indep        integer         independent combinations
  !>   rss              real(rk)        d . S . n by 3x3 algebra [MPa]
  !>   s_mat            real(rk)(3, 3)  the vertex as a matrix [MPa]
  !>   active5          real(rk)(5, 5)  Gram matrix of the 5-vectors of the systems at their
  !>                                    strength
  !>   lu, det, perm, status            its factorization and determinant
  !>   inside, spans, distinct  logical  the checks over all vertices
  !>   i, j, b          integer         vertices, system
  subroutine check_properties(name)

    character(len=*), intent(in) :: name

    integer :: n, num_indep, i, j, b, perm(5), status
    real(rk) :: strength(crys%num_slip), rss, s_mat(3, 3), active5(5, 5), lu(5, 5), det, p5(5)
    real(rk), allocatable :: vertices(:, :), scaled(:, :), reversed(:, :)
    logical :: inside, spans, distinct, deviatoric

    n = crys%num_slip
    call random_number(strength)
    strength = 100.0_rk + 100.0_rk*strength
    call yield_vertices_compute(crys%schmid_sym_c, strength, vertices, num_indep)
    call assert_true(size(vertices, 2) > 0, name // ': vertices found')

    inside = .true.
    spans = .true.
    deviatoric = .true.
    do i = 1, size(vertices, 2)
      ! Deviatoric; inside every hyperplane pair by explicit 3x3 algebra; the systems at their
      ! strength span the deviatoric space.
      deviatoric = deviatoric .and. abs(mandel_trace(vertices(:, i))) <= 1.0e-12_rk*maxval(strength)
      s_mat = matrix_from_mandel(vertices(:, i))
      active5 = 0.0_rk
      do b = 1, n
        rss = dot_product(crys%slip_dir_c(:, b), matmul(s_mat, crys%slip_normal_c(:, b)))
        if (abs(rss) > strength(b)*(1.0_rk + 1.0e-10_rk)) inside = .false.
        if (abs(abs(rss) - strength(b)) <= 1.0e-9_rk*strength(b)) then
          p5 = dev5_from_mandel(crys%schmid_sym_c(:, b))
          do j = 1, 5
            active5(:, j) = active5(:, j) + p5*p5(j)
          end do
        end if
      end do
      call lu_factor(active5, lu, perm, det, status)
      if (status /= 0 .or. .not. (abs(det) > 1.0e-12_rk*maxval(abs(active5))**5)) &
        & spans = .false.
    end do
    call assert_true(deviatoric, name // ': vertices deviatoric (Mandel trace 0)')
    call assert_true(inside, name // ': every vertex inside every hyperplane pair (3x3)')
    call assert_true(spans, name // ': the active systems of every vertex span 5 dimensions')

    ! No duplicates, and no vertex with its negative.
    distinct = .true.
    do i = 1, size(vertices, 2)
      do j = i + 1, size(vertices, 2)
        if (norm2(vertices(:, i) - vertices(:, j)) <= yield_tol_dup*minval(strength) .or. &
          & norm2(vertices(:, i) + vertices(:, j)) <= yield_tol_dup*minval(strength)) &
          & distinct = .false.
      end do
    end do
    call assert_true(distinct, name // ': no duplicate vertices or +- pairs')

    ! Scaling of the strengths scales the vertices (same order).
    call yield_vertices_compute(crys%schmid_sym_c, 3.0_rk*strength, scaled, num_indep)
    call assert_true(size(scaled, 2) == size(vertices, 2), name // ': scaled, same number')
    if (size(scaled, 2) == size(vertices, 2)) call assert_true( &
      & maxval(abs(scaled - 3.0_rk*vertices)) <= 1.0e-12_rk*3.0_rk*maxval(abs(vertices)), &
      & name // ': strengths x 3 give vertices x 3')

    ! Reversed order of the systems: the same set.
    call yield_vertices_compute(crys%schmid_sym_c(:, n:1:-1), strength(n:1:-1), reversed, &
      & num_indep)
    call assert_equal(size(reversed, 2), size(vertices, 2), name // ': reversed, same number')
    call assert_true(same_set(vertices, reversed, 1.0e-10_rk*maxval(strength)), &
      & name // ': reversed systems, the same set')

  end subroutine check_properties

  !> hcp basal + prismatic only spans 4 deviatoric dimensions: no vertex (header).
  !>
  !> Locals:
  !>   vertices   real(rk)(:, :)  computed vertices
  !>   num_indep  integer         independent combinations
  subroutine check_rank_deficient()

    real(rk), allocatable :: vertices(:, :)
    integer :: num_indep

    call crystal_init(crys, crystal_hcp, c_over_a=1.587_rk, families=[1, 2])
    call yield_vertices_compute(crys%schmid_sym_c, spread(100.0_rk, 1, 6), vertices, num_indep)
    call assert_equal(num_indep, 0, 'hcp basal + prismatic: no independent combination')
    call assert_equal(size(vertices, 2), 0, 'hcp basal + prismatic: no vertex')

  end subroutine check_rank_deficient

  !> Mandel 6-vectors of KD25 5-vectors (addendum_5vectors.md (F4)).
  !>
  !> Arguments:
  !>   ref5  (in)  real(rk)(5, :)  5-vectors
  !>   vec   (out, result) real(rk)(6, :)  Mandel 6-vectors
  !> Locals:
  !>   i     integer  column
  function mandel_from_ref(ref5) result(vec)

    real(rk), intent(in) :: ref5(:, :)
    real(rk) :: vec(6, size(ref5, 2))

    integer :: i

    do i = 1, size(ref5, 2)
      vec(1, i) = ref5(1, i)/sqrt(2.0_rk) - ref5(2, i)/sqrt(6.0_rk)
      vec(2, i) = -ref5(1, i)/sqrt(2.0_rk) - ref5(2, i)/sqrt(6.0_rk)
      vec(3, i) = 2.0_rk*ref5(2, i)/sqrt(6.0_rk)
      vec(4:6, i) = ref5(3:5, i)
    end do

  end function mandel_from_ref

  !> Whether every vector of a lies within tol of + or - a vector of b, and every vector of b
  !> within tol of + or - a vector of a (sets of +- pairs, either sign stored, any order).
  !>
  !> Arguments:
  !>   a, b  (in)  real(rk)(6, :)  the two sets
  !>   tol   (in)  real(rk)        distance tolerance
  !>   same  (out, result) logical
  !> Locals:
  !>   i     integer  vector
  function same_set(a, b, tol) result(same)

    real(rk), intent(in) :: a(:, :)
    real(rk), intent(in) :: b(:, :)
    real(rk), intent(in) :: tol
    logical :: same

    integer :: i

    same = .true.
    do i = 1, size(a, 2)
      same = same .and. min(distance_to_set(a(:, i), b), distance_to_set(-a(:, i), b)) <= tol
    end do
    do i = 1, size(b, 2)
      same = same .and. min(distance_to_set(b(:, i), a), distance_to_set(-b(:, i), a)) <= tol
    end do

  end function same_set

  !> Distance from a vector to the nearest column of a set.
  !>
  !> Arguments:
  !>   v     (in)  real(rk)(6)     vector
  !>   set   (in)  real(rk)(6, :)  set
  !>   dist  (out, result) real(rk)
  !> Locals:
  !>   j     integer  column
  function distance_to_set(v, set) result(dist)

    real(rk), intent(in) :: v(6)
    real(rk), intent(in) :: set(:, :)
    real(rk) :: dist

    integer :: j

    dist = huge(1.0_rk)
    do j = 1, size(set, 2)
      dist = min(dist, norm2(set(:, j) - v))
    end do

  end function distance_to_set

end program test_yield_vertices
