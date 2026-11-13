! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Periodic boundary conditions (ROADMAP step 8.3b; D71-D73): the drive pairs, the domain lengths
!> and the replicated constraint table v = T v_hat + g of the periodicity relations, the
!> consistent nodal weights of the zero-mean fluctuation, and the macroscopic averages.
!>
!> Theory : dev/theory/addendum_periodic_mpc.md §1.2, §3.3 (approved, revision 3): every
!>          relation (s, p, P) of $Periodicity (x_s = x_p + P o l) imposes, for each velocity
!>          component k, v_k(s) - v_k(p) = sum_j P_j l_j Fbar-dot_kj (R-P4 Eq. 73; Miehe & Koch
!>          2002, Eqs. 14, 17, 56-58, R-P2), with Fbar-dot the rate of the macroscopic
!>          deformation gradient on the reference periodicity vector and l the initial lengths;
!>          the macroscopic Cauchy stress is the volume average (R-P2); Fbar the volume average
!>          of F (R-P4 Eq. 68) and the drive pairs' separations (R-P4 Eq. 72).
!> Source : port of v2-fixed's periodicity (periodicity_mod1.f90::set_periodicity,
!>          periodicity_mod2.f90, periodicity_mod3.f90; conjugate_gradient_mod2.f90::prim_to_all
!>          "PBC"): the drive pair of axis j is the first relation of the file whose vector is
!>          sigma_j e_j (fill_pv_table, R-P4's reference node pair, Eq. 77); l_j = |x_j(Sd_j) -
!>          x_j(Pd_j)| (derive_gage_length); a free component Fbar-dot_kj is carried by the drive
!>          pair: Fbar-dot_kj = sigma_j (v_k(Sd_j) - v_k(Pd_j)) / l_j, v_k(Sd_j) freed
!>          (free_non_imposed_drive_nodes), each axis term with its own sign sigma_j (invert_drive
!>          after v2-fixed C70); diagonal vectors split into the axis terms
!>          (split_diagonal_periodicity_vectors; R-P4 Eq. 78); an imposed component gives the
!>          offset P_j l_j Fbar-dot_kj (apply_loading_as_offset); small strain (v2's strainrate,
!>          not in R-P4): Fbar-dot symmetric, one unknown per free shear pair {a, a+1} (cyclic),
!>          v_a(Sd_(a+1)), the other component of the pair scaled by l_j/l_k (C71, C72:
!>          link_symmetric_shear_components before the freeing, with the length factor).
!>          Deviations (addendum §3.3): the mesh must be periodic along all three axes (v2:
!>          undefined drive pairs otherwise, P8); chains (a primary that is a secondary, a drive
!>          node that is constrained) are resolved by substitution, so that no column of T is a
!>          constrained dof (Neper's meshes have none: then T is v2's); the rigid translation is
!>          fixed by the zero-mean fluctuation (Q-P3, constraint_mod::constraint_zero_mean).
!>
!> Policy (layout of mesh_mod and constraint_mod):
!>   - Components of Fbar and Fbar-dot are indexed c = 3 (i - 1) + j for ij.
!>   - The table is built identically on every rank (periodic_table is pure); the lengths come
!>     from the drive nodes' coordinates, summed over the ranks from their owners.
!>   - Collective routines: periodic_build, periodic_stress_average, periodic_fbar_average.
!>
!> Public routines:
!>   periodic_drive_pairs     the drive pair of each axis and its sign (pure)
!>   periodic_table           the replicated constraint table (pure; serial tests)
!>   periodic_build           the constraints of a periodic run (collective)
!>   periodic_stress_average  volume average of the Cauchy stress (collective)
!>   periodic_fbar_average    volume average of the deformation gradient (collective)
module periodic_mod

  use kinds_mod, only: rk, ik
  use constants_mod, only: sqrt2
  use error_mod, only: fatal_error
  use log_mod, only: log_info, log_warning
  use mesh_mod, only: mesh_t, node_loc_from_glob, find_sorted
  use element_mod, only: element_t, element_geometry_qpt
  use parallel_mod, only: par_sum
  use halo_mod, only: halo_forward, halo_reverse_add
  use constraint_mod, only: constraint_t, constraint_localize, sort_index, unique_sorted

  implicit none

  private

  public :: periodic_drive_pairs, periodic_table, periodic_build, periodic_stress_average, &
    & periodic_fbar_average

  ! Component names of Fbar (index 3 (i - 1) + j), for messages.
  character(len=2), parameter :: comp_name(9) = ['11', '12', '13', '21', '22', '23', '31', &
    & '32', '33']
  ! Axis letters, for messages.
  character(len=1), parameter :: axis_letter(3) = ['x', 'y', 'z']

contains

  !> The drive pair of each axis j: the first relation of the file whose vector is sigma_j e_j,
  !> sigma_j = +-1 (v2 fill_pv_table). A mesh without such a relation for some axis is not
  !> periodic along it: an error (addendum §3.3).
  !>
  !> Arguments:
  !>   per_secondary, per_primary  (in)  integer(ik)(:)     relations' node ids
  !>   per_vect    (in)   integer(ik)(:, :)  (3, :) their vectors P
  !>   drive_node  (out)  integer(ik)(2, 3)  secondary (1) and primary (2) node of each axis
  !>   sigma       (out)  real(rk)(3)        sign of each drive vector
  !> Locals:
  !>   i, j        integer(ik)  relation, axis
  subroutine periodic_drive_pairs(per_secondary, per_primary, per_vect, drive_node, sigma)

    integer(ik), intent(in) :: per_secondary(:), per_primary(:), per_vect(:, :)
    integer(ik), intent(out) :: drive_node(2, 3)
    real(rk), intent(out) :: sigma(3)

    integer(ik) :: i, j

    drive_node = 0_ik
    sigma = 0.0_rk
    do i = 1_ik, int(size(per_secondary), ik)
      if (sum(abs(per_vect(:, i))) /= 1_ik) cycle
      j = findloc(abs(per_vect(:, i)), 1_ik, dim=1)
      if (drive_node(1, j) > 0_ik) cycle
      drive_node(:, j) = [per_secondary(i), per_primary(i)]
      sigma(j) = real(per_vect(j, i), rk)
    end do
    do j = 1_ik, 3_ik
      if (drive_node(1, j) == 0_ik) call fatal_error('periodic', 'the mesh is not periodic ' // &
        & 'along ' // axis_letter(j) // ' ($Periodicity has no relation with vector +-e_' &
        & // axis_letter(j) // '): periodic conditions need a mesh periodic along all ' // &
        & 'three axes (neper -T -periodicity all)')
    end do

  end subroutine periodic_drive_pairs

  !> The replicated constraint table of the periodicity relations (module header): one row per
  !> secondary dof but the freed drive dofs, with the primary, the drive columns of the free
  !> components and the rate offsets of the imposed ones; chains resolved by substitution;
  !> repeated columns merged. Pure in effect (errors are fatal and the same on every rank).
  !>
  !> Arguments:
  !>   per_secondary, per_primary  (in)  integer(ik)(:)  relations' node ids
  !>   per_vect    (in)  integer(ik)(:, :)  (3, :) their vectors P
  !>   drive_node  (in)  integer(ik)(2, 3)  drive pairs (periodic_drive_pairs)
  !>   sigma       (in)  real(rk)(3)        their signs
  !>   lengths     (in)  real(rk)(3)        domain lengths l_j [mm]
  !>   small       (in)  logical            small strain (symmetric Fbar-dot)
  !>   imposed     (in)  logical(9)         imposed components of Fbar-dot (small strain: ij and
  !>                                        ji alike)
  !>   cons        (inout) type(constraint_t)  glob_* and all_cols set
  !> Locals:
  !>   max_cols    integer(ik)        columns of a raw row (1 + 2 per axis term)
  !>   n_raw       integer(ik)        raw rows
  !>   raw_dof, raw_prim  integer(ik)(:)  constrained dof and primary dof of each raw row
  !>   raw_ncol    integer(ik)(:)     columns of each raw row
  !>   raw_col     integer(ik)(:, :)  (max_cols, n_raw) column dofs
  !>   raw_coef    real(rk)(:, :)     their coefficients
  !>   raw_rate    real(rk)(:, :)     (9, n_raw) offset per unit rate [mm]
  !>   freed       integer(ik)(9)     freed drive dofs (0: none)
  !>   order       integer(ik)(:)     raw rows by dof
  !>   s_dof       integer(ik)(:)     row dofs, increasing
  !>   i, k, j, c  integer(ik)        relation, component, axis, rate index
  !>   kk, jj      integer(ik)        component and axis of the drive dof carrying a component
  !>   n, r, q, e  integer(ik)        counters, rows, columns
  !>   f           real(rk)           coefficient of a drive term
  !>   w_col, w_coef  (:)             columns and coefficients of a row being resolved
  !>   w_rate      real(rk)(9)        its rate offsets
  !>   pos, steps  integer(ik)        a row found among the rows, substitutions made
  !>   prim        integer(ik)        top primary dof
  !>   srt         integer(ik)(:)     columns sorted
  !>   m_col, m_coef  (:)             merged columns and coefficients
  !>   n_tot       integer(ik)        columns of the table
  subroutine periodic_table(per_secondary, per_primary, per_vect, drive_node, sigma, lengths, &
    & small, imposed, cons)

    integer(ik), intent(in) :: per_secondary(:), per_primary(:), per_vect(:, :)
    integer(ik), intent(in) :: drive_node(2, 3)
    real(rk), intent(in) :: sigma(3), lengths(3)
    logical, intent(in) :: small, imposed(9)
    type(constraint_t), intent(inout) :: cons

    integer(ik), parameter :: max_cols = 7_ik
    integer(ik) :: n_raw, freed(9), i, k, j, c, kk, jj, n, r, q, e, pos, steps, prim, n_tot
    integer(ik), allocatable :: raw_dof(:), raw_prim(:), raw_ncol(:), raw_col(:, :), order(:), &
      & s_dof(:), w_col(:), srt(:), m_col(:)
    real(rk), allocatable :: raw_coef(:, :), raw_rate(:, :), w_coef(:), m_coef(:)
    real(rk) :: f, w_rate(9)
    type :: row_t
      integer(ik), allocatable :: col(:)
      real(rk), allocatable :: coef(:)
    end type row_t
    type(row_t), allocatable :: rows(:)

    ! Freed drive dofs: the dof carrying each free component (v2 free_non_imposed_drive_nodes).
    freed = 0_ik
    do k = 1_ik, 3_ik
      do j = 1_ik, 3_ik
        c = 3_ik*(k - 1_ik) + j
        if (imposed(c)) cycle
        call carrier(k, j, kk, jj)
        freed(c) = 3_ik*(drive_node(1, jj) - 1_ik) + kk
      end do
    end do

    ! Raw rows: every secondary dof but the freed ones, v_k(s) = v_k(p) + sum_j P_j l_j
    ! Fbar-dot_kj, each free Fbar-dot_kj by its drive pair, each imposed one as an offset.
    n_raw = 3_ik*int(size(per_secondary), ik)
    allocate (raw_dof(n_raw), raw_prim(n_raw), raw_ncol(n_raw), raw_col(max_cols, n_raw), &
      & raw_coef(max_cols, n_raw), raw_rate(9, n_raw))
    n = 0_ik
    do i = 1_ik, int(size(per_secondary), ik)
      do k = 1_ik, 3_ik
        if (any(freed == 3_ik*(per_secondary(i) - 1_ik) + k)) cycle
        n = n + 1_ik
        raw_dof(n) = 3_ik*(per_secondary(i) - 1_ik) + k
        raw_prim(n) = 3_ik*(per_primary(i) - 1_ik) + k
        raw_ncol(n) = 1_ik
        raw_col(1, n) = raw_prim(n)
        raw_coef(1, n) = 1.0_rk
        raw_rate(:, n) = 0.0_rk
        do j = 1_ik, 3_ik
          if (per_vect(j, i) == 0_ik) cycle
          c = 3_ik*(k - 1_ik) + j
          if (imposed(c)) then
            raw_rate(c, n) = raw_rate(c, n) + real(per_vect(j, i), rk)*lengths(j)
          else
            ! The drive pair of axis jj, component kk, scaled by sigma_jj / l_jj (and l_j / l_k
            ! for the relinked shear slot of small strain, C72).
            call carrier(k, j, kk, jj)
            f = real(per_vect(j, i), rk)*lengths(j)*sigma(jj)/lengths(jj)
            raw_col(raw_ncol(n) + 1_ik:raw_ncol(n) + 2_ik, n) = &
              & [3_ik*(drive_node(1, jj) - 1_ik) + kk, 3_ik*(drive_node(2, jj) - 1_ik) + kk]
            raw_coef(raw_ncol(n) + 1_ik:raw_ncol(n) + 2_ik, n) = [f, -f]
            raw_ncol(n) = raw_ncol(n) + 2_ik
          end if
        end do
      end do
    end do

    ! Rows by dof; a node that is the secondary of two relations is an error.
    order = sort_index(raw_dof(1:n))
    s_dof = raw_dof(order)
    do r = 2_ik, n
      if (s_dof(r) == s_dof(r - 1_ik)) call fatal_error('periodic_table', 'node ' // &
        & int_text((s_dof(r) - 1_ik)/3_ik + 1_ik) // ' is the secondary of two periodicity ' // &
        & 'relations')
    end do

    ! Each row resolved: a column that is itself a constrained dof is replaced by its row
    ! (coefficients and rate offsets scaled); the top primary followed likewise.
    allocate (rows(n), cons%glob_rate_coef(9, n), cons%glob_prim(n))
    do r = 1_ik, n
      q = order(r)
      w_col = raw_col(1:raw_ncol(q), q)
      w_coef = raw_coef(1:raw_ncol(q), q)
      w_rate = raw_rate(:, q)
      steps = 0_ik
      e = 1_ik
      do while (e <= int(size(w_col), ik))
        pos = find_sorted(s_dof, w_col(e))
        if (pos == 0_ik) then
          e = e + 1_ik
          cycle
        end if
        steps = steps + 1_ik
        if (steps > n) call fatal_error('periodic_table', 'the periodicity relations form ' // &
          & 'a cycle (node ' // int_text((s_dof(r) - 1_ik)/3_ik + 1_ik) // ')')
        f = w_coef(e)
        q = order(pos)
        w_rate = w_rate + f*raw_rate(:, q)
        w_col = [w_col(1:e - 1_ik), w_col(e + 1_ik:), raw_col(1:raw_ncol(q), q)]
        w_coef = [w_coef(1:e - 1_ik), w_coef(e + 1_ik:), f*raw_coef(1:raw_ncol(q), q)]
      end do
      prim = raw_prim(order(r))
      steps = 0_ik
      do
        pos = find_sorted(s_dof, prim)
        if (pos == 0_ik) exit
        steps = steps + 1_ik
        if (steps > n) call fatal_error('periodic_table', 'the periodicity relations form ' // &
          & 'a cycle (node ' // int_text((s_dof(r) - 1_ik)/3_ik + 1_ik) // ')')
        prim = raw_prim(order(pos))
      end do
      ! Repeated columns merged (coefficients summed; a zero sum is kept, as the column may be
      ! the top primary).
      srt = sort_index(w_col)
      allocate (m_col(size(w_col)), m_coef(size(w_col)))
      q = 0_ik
      do e = 1_ik, int(size(w_col), ik)
        if (q > 0_ik) then
          if (m_col(q) == w_col(srt(e))) then
            m_coef(q) = m_coef(q) + w_coef(srt(e))
            cycle
          end if
        end if
        q = q + 1_ik
        m_col(q) = w_col(srt(e))
        m_coef(q) = w_coef(srt(e))
      end do
      rows(r)%col = m_col(1:q)
      rows(r)%coef = m_coef(1:q)
      deallocate (m_col, m_coef)
      cons%glob_rate_coef(:, r) = w_rate
      cons%glob_prim(r) = prim
    end do

    ! The table in the layout of constraint_t (rows in increasing dof, constant offsets 0).
    cons%num_rows_glob = n
    n_tot = 0_ik
    do r = 1_ik, n
      n_tot = n_tot + int(size(rows(r)%col), ik)
    end do
    allocate (cons%glob_dof(n), cons%glob_start(n + 1_ik), cons%glob_col(n_tot), &
      & cons%glob_coef(n_tot), cons%glob_offset(n))
    cons%glob_dof = s_dof
    cons%glob_offset = 0.0_rk
    e = 0_ik
    do r = 1_ik, n
      cons%glob_start(r) = e + 1_ik
      cons%glob_col(e + 1_ik:e + int(size(rows(r)%col), ik)) = rows(r)%col
      cons%glob_coef(e + 1_ik:e + int(size(rows(r)%col), ik)) = rows(r)%coef
      e = e + int(size(rows(r)%col), ik)
    end do
    cons%glob_start(n + 1_ik) = e + 1_ik
    cons%all_cols = unique_sorted(cons%glob_col)

  contains

    !> The drive dof that carries Fbar-dot_kj: component kk of the drive pair of axis jj
    !> (deformation_gradient_rate: kk = k, jj = j; small strain: the diagonal and the pair
    !> {a, a+1} cyclic by component a of axis a + 1, v2's relinking, C71, C72).
    !>
    !> Arguments:
    !>   k, j    (in)   integer(ik)  the component kj
    !>   kk, jj  (out)  integer(ik)  component and axis of its drive dof
    subroutine carrier(k, j, kk, jj)

      integer(ik), intent(in) :: k, j
      integer(ik), intent(out) :: kk, jj

      kk = k
      jj = j
      if (small .and. j /= k .and. j /= mod(k, 3_ik) + 1_ik) then
        kk = j
        jj = k
      end if

    end subroutine carrier

  end subroutine periodic_table

  !> The constraints of a periodic run (collective): the drive pairs, the lengths l_j from the
  !> drive nodes' initial coordinates (v2 derive_gage_length), the replicated table, its local
  !> rows, the preconditioner choice, the consistent nodal weights w_a = int N_a dOmega_0 of the
  !> initial configuration, V_0 and W_j = sum_a w_a X_a,j; a warning for a free rigid rotation
  !> rate (both ij and ji free with deformation_gradient_rate, Q-P2).
  !>
  !> Arguments:
  !>   cons     (out)  type(constraint_t)  the constraints
  !>   mesh     (in)   type(mesh_t)        rank-local mesh (relations, initial coordinates)
  !>   elt      (in)   type(element_t)     element descriptor
  !>   small    (in)   logical             small strain kinematics
  !>   imposed  (in)   logical(9)          imposed components of Fbar-dot
  !>   exact_diagonal  (in)  logical       the exact diagonal of T^T K T (else the
  !>                                       lumped diagonal)
  !> Locals:
  !>   sigma    real(rk)(3)        drive signs
  !>   part, whole  real(rk)(3, 2, 3)  drive nodes' coordinates: owned, over the ranks [mm]
  !>   flat     real(rk)(18)       the sums, flattened
  !>   lengths  real(rk)(3)        l_j [mm]
  !>   j, end_pair, a, i, k  integer  axis, pair end, element node, components
  !>   i_loc    integer(ik)        local node
  !>   i_elt, i_qpt  integer(ik)   element, point
  !>   dshape_s real(rk)(:, :)     shape derivatives (unused) [1/mm]
  !>   det_jac  real(rk)           det J [mm^3]
  !>   sums, sums_all  real(rk)(4)  V_0 and W over the owned nodes, over the ranks
  !>   num_cols integer(ik)        columns of the table
  !>   text     character(240)     message
  subroutine periodic_build(cons, mesh, elt, small, imposed, exact_diagonal)

    type(constraint_t), intent(out) :: cons
    type(mesh_t), intent(in) :: mesh
    type(element_t), intent(in) :: elt
    logical, intent(in) :: small, imposed(9), exact_diagonal

    real(rk) :: sigma(3), part(3, 2, 3), whole(3, 2, 3), flat(18), lengths(3), det_jac
    real(rk) :: sums(4), sums_all(4)
    real(rk), allocatable :: dshape_s(:, :)
    integer :: j, end_pair, a, i, k
    integer(ik) :: i_loc, i_elt, i_qpt, num_cols
    character(len=240) :: text

    if (mesh%num_periodicity == 0_ik) call fatal_error('periodic', 'boundary_conditions ' // &
      & 'periodic needs a mesh with periodicity relations ($Periodicity; neper -T ' // &
      & '-periodicity all)')

    ! Drive pairs and their nodes' coordinates (from the owners).
    call periodic_drive_pairs(mesh%per_secondary, mesh%per_primary, mesh%per_vect, &
      & cons%drive_node, sigma)
    part = 0.0_rk
    do j = 1, 3
      do end_pair = 1, 2
        i_loc = node_loc_from_glob(mesh, cons%drive_node(end_pair, j))
        cons%drive_loc(end_pair, j) = 0_ik
        if (i_loc > 0_ik .and. i_loc <= mesh%num_nodes_owned) then
          cons%drive_loc(end_pair, j) = i_loc
          part(:, end_pair, j) = mesh%coo_s(:, i_loc)
        end if
      end do
    end do
    call par_sum(reshape(part, [18]), flat)
    whole = reshape(flat, [3, 2, 3])

    ! Lengths: the drive pairs' separations along their axes (v2 derive_gage_length).
    do j = 1, 3
      lengths(j) = abs(whole(j, 1, j) - whole(j, 2, j))
      if (.not. lengths(j) > 0.0_rk) call fatal_error('periodic', 'the drive pair of axis ' // &
        & axis_letter(j) // ' has zero length')
      cons%drive_factor(j) = sigma(j)/lengths(j)
    end do
    write (text, '(a,3(es11.4,a),6(i0,a))') 'Periodic conditions: lengths ', lengths(1), ', ', &
      & lengths(2), ', ', lengths(3), ' mm; drive pairs (secondary, primary) x (', &
      & cons%drive_node(1, 1), ', ', cons%drive_node(2, 1), '), y (', cons%drive_node(1, 2), &
      & ', ', cons%drive_node(2, 2), '), z (', cons%drive_node(1, 3), ', ', &
      & cons%drive_node(2, 3), ')'
    call log_info(trim(text), 2)

    ! A free rigid rotation rate (Q-P2): allowed, with a warning.
    if (.not. small) then
      do i = 1, 3
        do k = i + 1, 3
          if (imposed(3*(i - 1) + k) .or. imposed(3*(k - 1) + i)) cycle
          call log_warning('Periodic conditions: Fbar-dot ' // comp_name(3*(i - 1) + k) // &
            & ' and ' // comp_name(3*(k - 1) + i) // ' are both free: the rigid rotation ' // &
            & 'rate in that plane is not fixed (impose one of them, or use small_strain).', 2)
        end do
      end do
    end if

    ! The table, its local rows and the preconditioner.
    call periodic_table(mesh%per_secondary, mesh%per_primary, mesh%per_vect, cons%drive_node, &
      & sigma, lengths, small, imposed, cons)
    cons%active = .true.
    cons%diag_exact = exact_diagonal
    call constraint_localize(cons, mesh)
    num_cols = int(size(cons%all_cols), ik)
    call log_info('Periodic conditions: ' // int_text(cons%num_rows_glob) // &
      & ' constrained dofs, ' // int_text(num_cols) // ' independent dofs referenced', 2)

    ! Consistent nodal weights of the initial configuration, complete at the owned and ghost
    ! nodes; V_0 and W.
    allocate (cons%node_weight(mesh%num_nodes_loc), dshape_s(3, elt%num_nodes))
    cons%node_weight = 0.0_rk
    do i_elt = 1_ik, mesh%num_elts
      do i_qpt = 1_ik, elt%num_qpts
        call element_geometry_qpt(elt, mesh%coo_s(:, mesh%elt_nodes_loc(:, i_elt)), i_qpt, &
          & dshape_s, det_jac)
        do a = 1, int(elt%num_nodes)
          i_loc = mesh%elt_nodes_loc(a, i_elt)
          cons%node_weight(i_loc) = cons%node_weight(i_loc) + &
            & elt%qpt_weight(i_qpt)*det_jac*elt%shape_qpt(a, i_qpt)
        end do
      end do
    end do
    call halo_reverse_add(mesh, cons%node_weight, 1_ik)
    call halo_forward(mesh, cons%node_weight, 1_ik)
    sums = 0.0_rk
    do i_loc = 1_ik, mesh%num_nodes_owned
      sums = sums + cons%node_weight(i_loc)*[1.0_rk, mesh%coo_s(:, i_loc)]
    end do
    call par_sum(sums, sums_all)
    cons%volume0 = sums_all(1)
    cons%moment0 = sums_all(2:4)
    cons%zero_mean = .true.

  end subroutine periodic_build

  !> Volume average of the Cauchy stress over the current configuration (R-P2; addendum §3.3,
  !> Q-P11): sum over the elements and points of w_q det J sigma_q over the volume. Collective.
  !>
  !> Arguments:
  !>   mesh      (in)  type(mesh_t)       rank-local mesh
  !>   elt       (in)  type(element_t)    element descriptor
  !>   coo       (in)  real(rk)(:, :)     (3, num_nodes_loc) current coordinates [mm]
  !>   stress_s  (in)  real(rk)(:, :, :)  (6, num_qpts, num_elts) Cauchy stress, Mandel, sample
  !>                                      frame [MPa]
  !>   avg       (out, result) real(rk)(6)  the average, tensor components 11 22 33 23 13 12
  !>                                      [MPa]
  !> Locals:
  !>   part, whole  real(rk)(7)    sum of w det J sigma (Mandel) and the volume: here, all ranks
  !>   dshape_s     real(rk)(:, :) shape derivatives (unused) [1/mm]
  !>   det_jac      real(rk)       det J [mm^3]
  !>   i_elt, i_qpt integer(ik)    element, point
  function periodic_stress_average(mesh, elt, coo, stress_s) result(avg)

    type(mesh_t), intent(in) :: mesh
    type(element_t), intent(in) :: elt
    real(rk), intent(in) :: coo(:, :)
    real(rk), intent(in) :: stress_s(:, :, :)
    real(rk) :: avg(6)

    real(rk) :: part(7), whole(7), det_jac
    real(rk), allocatable :: dshape_s(:, :)
    integer(ik) :: i_elt, i_qpt

    allocate (dshape_s(3, elt%num_nodes))
    part = 0.0_rk
    do i_elt = 1_ik, mesh%num_elts
      do i_qpt = 1_ik, elt%num_qpts
        call element_geometry_qpt(elt, coo(:, mesh%elt_nodes_loc(:, i_elt)), i_qpt, dshape_s, &
          & det_jac)
        part = part + elt%qpt_weight(i_qpt)*det_jac*[stress_s(:, i_qpt, i_elt), 1.0_rk]
      end do
    end do
    call par_sum(part, whole)

    ! Mandel to tensor components.
    avg = whole(1:6)/whole(7)
    avg(4:6) = avg(4:6)/sqrt2

  end function periodic_stress_average

  !> Volume average of the deformation gradient over the initial configuration (R-P4 Eq. 68):
  !> (1/V_0) sum over the elements and points of w_q det J_0 dx/dX. Collective.
  !>
  !> Arguments:
  !>   mesh   (in)  type(mesh_t)     rank-local mesh (initial coordinates coo_s)
  !>   elt    (in)  type(element_t)  element descriptor
  !>   coo    (in)  real(rk)(:, :)   (3, num_nodes_loc) current coordinates [mm]
  !>   f_avg  (out, result) real(rk)(3, 3)  the average F_ij [-]
  !> Locals:
  !>   part, whole  real(rk)(10)   sum of w det J_0 F (column major) and the volume
  !>   dshape_s     real(rk)(:, :) dN_a/dX_j, initial configuration [1/mm]
  !>   det_jac      real(rk)       det J_0 [mm^3]
  !>   f_pt         real(rk)(3, 3) F at a point [-]
  !>   i_elt, i_qpt integer(ik)    element, point
  function periodic_fbar_average(mesh, elt, coo) result(f_avg)

    type(mesh_t), intent(in) :: mesh
    type(element_t), intent(in) :: elt
    real(rk), intent(in) :: coo(:, :)
    real(rk) :: f_avg(3, 3)

    real(rk) :: part(10), whole(10), det_jac, f_pt(3, 3)
    real(rk), allocatable :: dshape_s(:, :)
    integer(ik) :: i_elt, i_qpt

    allocate (dshape_s(3, elt%num_nodes))
    part = 0.0_rk
    do i_elt = 1_ik, mesh%num_elts
      do i_qpt = 1_ik, elt%num_qpts
        call element_geometry_qpt(elt, mesh%coo_s(:, mesh%elt_nodes_loc(:, i_elt)), i_qpt, &
          & dshape_s, det_jac)
        f_pt = matmul(coo(:, mesh%elt_nodes_loc(:, i_elt)), transpose(dshape_s))
        part = part + elt%qpt_weight(i_qpt)*det_jac*[reshape(f_pt, [9]), 1.0_rk]
      end do
    end do
    call par_sum(part, whole)
    f_avg = reshape(whole(1:9), [3, 3])/whole(10)

  end function periodic_fbar_average

  !> An integer as text, without blanks.
  !>
  !> Arguments:
  !>   value  (in)  integer(ik)  the integer
  !>   text   (out, result) character(:), allocatable
  !> Locals:
  !>   buffer  character(24)  formatted value
  pure function int_text(value) result(text)

    integer(ik), intent(in) :: value
    character(len=:), allocatable :: text

    character(len=24) :: buffer

    write (buffer, '(i0)') value
    text = trim(buffer)

  end function int_text

end module periodic_mod
