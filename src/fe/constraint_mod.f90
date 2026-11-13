! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Linear constraints between nodal velocity dofs, v = T v_hat + g (ROADMAP step 8.3a; D71, D73,
!> D74): multi-point constraints (MPCs, step 8.3a) and periodicity (step 8.3b, whose table is
!> built by periodic_mod).
!>
!> Theory : dev/theory/addendum_periodic_mpc.md §3.1-§3.2 (approved, revision 3): master-slave
!>          elimination, the constraint u = T u_hat + g and the reduced system T^T K T u_hat =
!>          T^T (f - K g) (Felippa, Introduction to Finite Element Methods, Ch. 8, §8.3,
!>          Eqs. 8.26-8.27 and 8.29-8.32, R-P1; Villette et al., in preparation, Eqs. 80-89 and
!>          Algorithm 1, R-P4).
!> Source : port of v2-fixed's MPCs: multi_point_constraints_mod1.f90::read_mpcs (each node of a
!>          node set but the first is a secondary of the first: v_s = coeff v_p + offset; a zero
!>          coefficient is an error), multi_point_constraints_mod2.f90::calc_mpcs (chains
!>          resolved to the top primary, coefficients multiplied and offsets composed) and
!>          ::propagate_bcs_mpcs (velocity conditions carried to the group, after CHANGES.md C66:
!>          the values given in a group averaged, a conflict larger than 1e-10 of the largest
!>          condition an error), conjugate_gradient_mod2.f90::prim_to_all / all_to_prim (T and
!>          T^T, with the coefficients after C67) and conjugate_gradient_mod1.f90::
!>          assemble_diagonals("MPC") (the diagonal of T^T K T, element by element, which since
!>          C67 is the exact one: addendum §1.1), driver_uniaxial_control_mod1.f90 l. 94-103 (the
!>          initial velocity projected onto the constraints). Deviations (addendum §3.2, D73): a
!>          dof that two lines make a secondary is an error (v2: the later line wins) unless the
!>          two relations are identical (then one is kept), an MPC
!>          cycle is an error (v2 would not stop), a node set that is not in the mesh is an error
!>          (v2 ignores the line), a node set with fewer than two nodes is a warning.
!>
!> Policy (layout of mesh_mod; parallel transfer: D74, addendum §3.1, constraint_exchange_mod):
!>   - The constraint table is replicated: every rank builds the same rows from replicated data
!>     (the MPC lines and mesh_t%nsets_glob). A row is a constrained global dof (3 (node id - 1)
!>     + k) with its columns (independent global dofs, never constrained ones: chains are
!>     resolved), their coefficients, and an offset [mm/s]. Rows are in increasing global dof.
!>   - Each rank keeps the rows of its local (owned and ghost) constrained dofs. A column is read
!>     from the rank's owned entry when it owns it, otherwise from its slot of the exchange
!>     (constraint_exchange_mod), whose set-up finds its owner.
!>   - A reduced vector (T^T applied, or a vector of independent dofs) is a full dof vector with
!>     zeros at the constrained dofs, so the distributed dot products over the owned dofs are
!>     unchanged (R-P4 §3.2.2).
!>   - T (constraint_expand): every local constrained dof, owned or ghost, is set from its
!>     columns; the unconstrained ghost dofs must be valid already (halo_forward), so that the
!>     element product needs no further exchange.
!>   - T^T (constraint_reduce): each owned constrained dof adds coef times its value to its
!>     columns (owned directly, others through the exchange's reverse add), then is set to zero.
!>   - Offsets (step 8.3b): g_s = offset0_s + sum_c a_sc rate_c, with rate the nine imposed
!>     rates of the macroscopic deformation gradient Fbar-dot of the step (index 3 (i - 1) + j
!>     for ij; constraint_set_rates); MPC offsets are constant (a = 0).
!>   - Preconditioner (step 8.3b, Q-P9, addendum §3.1): the exact diagonal of T^T K T, or the
!>     lumped periodic diagonal (v2's pre_cond 1: the element diagonals of every dof summed onto the top
!>     primary of its row, the "prim" column, the drive couplings left out; diag_exact false:
!>     constraint_diagonal_lumped, chosen by the caller).
!>   - Rigid translation of periodic conditions (step 8.3b, Q-P3, addendum §3.3): with the
!>     consistent nodal weights w_a = int N_a dOmega_0 (node_weight, owned and ghost), the
!>     weighted mean of each velocity component k not fixed by a velocity condition is
!>     sum_a w_a v_a,k = sum_j Fbar-dot_kj W_j, W_j = sum_a w_a X_a,j (constraint_zero_mean);
!>     Fbar-dot from the drive pairs (constraint_macro_rate, R-P4 Eq. 72 in rate form).
!>   - Collective routines: constraint_build_mpc, constraint_localize, constraint_expand,
!>     constraint_reduce, constraint_matvec, constraint_diagonal, constraint_propagate_bcs,
!>     constraint_project, constraint_macro_rate, constraint_zero_mean, constraint_check_bcs
!>     (every rank calls them, in the same order).
!>
!> Public types:
!>   mpc_spec_t    one MPC relation of the configuration (components at the type)
!>   constraint_t  the constraint table and its local rows (components at the type)
!> Public routines:
!>   constraint_table_mpc      the replicated table of MPC lines (pure; serial tests)
!>   constraint_build_mpc      the table, its local rows and the exchange (collective)
!>   constraint_localize       local rows of a replicated table and the exchange (collective)
!>   constraint_set_rates      offsets of the imposed macroscopic rates of a step
!>   constraint_expand         v = T v_hat (+ g) on the local constrained dofs
!>   constraint_reduce         y = T^T y on the independent dofs, zero at the constrained ones
!>   constraint_matvec         y = T^T K T x, zero at the constrained and prescribed dofs
!>   constraint_diagonal       exact diagonal of T^T K T
!>   constraint_diagonal_lumped    the lumped periodic diagonal (v2's pre_cond 1)
!>   constraint_propagate_bcs  velocity conditions carried through the MPC groups (C66)
!>   constraint_project        velocity projected onto the constraints, v = T v + g
!>   constraint_check_bcs      periodic: no velocity condition on a constrained dof; the
!>                             components fixed by conditions
!>   constraint_macro_rate     periodic: Fbar-dot from the drive pairs' velocities
!>   constraint_zero_mean      periodic: the velocity shifted to the zero-mean fluctuation
!>   sort_index, unique_sorted  integer sorting helpers (also used by periodic_mod)
module constraint_mod

  use kinds_mod, only: rk, ik
  use error_mod, only: fatal_error, internal_error
  use log_mod, only: log_info, log_warning
  use mesh_mod, only: mesh_t, node_loc_from_glob, find_sorted
  use mesh_file_mod, only: mesh_nset_t
  use parallel_mod, only: par_max, par_any, par_sum, par_log_flush
  use halo_mod, only: halo_forward, halo_reverse_add
  use ebe_mod, only: ebe_matvec
  use constraint_exchange_mod, only: constraint_exchange_t, constraint_exchange_setup, &
    & constraint_forward, constraint_reverse_add

  implicit none

  private

  public :: mpc_spec_t, constraint_t
  public :: constraint_table_mpc, constraint_build_mpc, constraint_expand, constraint_reduce, &
    & constraint_matvec, constraint_diagonal, constraint_propagate_bcs, constraint_project
  public :: constraint_localize, constraint_set_rates, constraint_diagonal_lumped, constraint_check_bcs, &
    & constraint_macro_rate, constraint_zero_mean
  public :: sort_index, unique_sorted

  !> Conflict tolerance of velocity conditions in an MPC group, relative to the largest
  !> condition (v2-fixed C66).
  real(rk), parameter :: conflict_rel_tol = 1.0e-10_rk

  ! Axis letters (messages).
  character(len=1), parameter :: axis_letter(3) = ['x', 'y', 'z']

  !> One MPC relation: every node s of the node set but the first (p) gets
  !> v_s,comp = coeff v_p,comp + offset.
  !>
  !> Components:
  !>   nset    character(:), allocatable  node-set label
  !>   comp    integer                    velocity component (1 x, 2 y, 3 z)
  !>   coeff   real(rk)                   coefficient (nonzero) [-]
  !>   offset  real(rk)                   offset [mm/s]
  !>   line    integer(ik)                line in the configuration file (messages)
  type :: mpc_spec_t
    character(len=:), allocatable :: nset
    integer :: comp = 0
    real(rk) :: coeff = 1.0_rk
    real(rk) :: offset = 0.0_rk
    integer(ik) :: line = 0_ik
  end type mpc_spec_t

  !> The constraint table (replicated) and this rank's local rows.
  !>
  !> Components:
  !>   active          logical          constraints are present
  !>   num_rows_glob   integer(ik)      rows of the whole table
  !>   glob_dof        integer(ik)(num_rows_glob)      constrained global dof of each row,
  !>                                    increasing
  !>   glob_start      integer(ik)(num_rows_glob + 1)  columns of row r: glob_start(r) ..
  !>                                    glob_start(r + 1) - 1
  !>   glob_col        integer(ik)(:)   independent global dof of each column
  !>   glob_coef       real(rk)(:)      coefficient of each column [-]
  !>   glob_offset     real(rk)(num_rows_glob)  constant offset of each row [mm/s]
  !>   glob_rate_coef  real(rk)(9, num_rows_glob)  offset of each row per unit imposed rate of
  !>                                    Fbar-dot (periodic; 0 for MPCs) [mm]
  !>   glob_prim       integer(ik)(num_rows_glob)  top primary dof of each row (one of its
  !>                                    columns; the lumped periodic diagonal)
  !>   all_cols        integer(ik)(:)   distinct columns of the table, increasing
  !>   num_rows        integer(ik)      local rows (constrained local dofs, owned or ghost)
  !>   row_dof         integer(ik)(num_rows)      local dof of each local row, increasing
  !>   row_start       integer(ik)(num_rows + 1)  columns of local row r
  !>   col_coef        real(rk)(:)      coefficient of each local column [-]
  !>   col_owned       integer(ik)(:)   owned local dof of the column, 0 if not owned here
  !>   col_slot        integer(ik)(:)   exchange slot of the column, 0 if owned here
  !>   col_local       integer(ik)(:)   local dof (owned or ghost) of the column, 0 if not local
  !>   col_glob        integer(ik)(:)   global dof of the column
  !>   row_offset      real(rk)(num_rows)  offset of each local row (with the rates set) [mm/s]
  !>   row_offset0     real(rk)(num_rows)  its constant part [mm/s]
  !>   row_rate_coef   real(rk)(9, num_rows)  its part per unit imposed rate [mm]
  !>   row_prim        integer(ik)(num_rows)  column entry of the row's top primary
  !>   row_of_dof      integer(ik)(3 num_nodes_loc)  local row of each local dof, 0 if the dof
  !>                                    is not constrained
  !>   is_constrained  logical(3 num_nodes_loc)      the dof is constrained
  !>   xch             type(constraint_exchange_t)   slots and lists of the transfer
  !>   diag_exact      logical          preconditioner: the exact diagonal of T^T K T (else
  !>                                    the lumped periodic diagonal)
  !>   zero_mean       logical          the rigid translation is fixed by the zero-mean
  !>                                    fluctuation (periodic)
  !>   mean_fixed      logical(3)       components fixed by velocity conditions (not shifted)
  !>   node_weight     real(rk)(num_nodes_loc)  consistent nodal weights int N_a dOmega_0, owned
  !>                                    and ghost [mm^3]
  !>   volume0         real(rk)         initial volume sum_a w_a [mm^3]
  !>   moment0         real(rk)(3)      W_j = sum_a w_a X_a,j, initial coordinates [mm^4]
  !>   drive_node      integer(ik)(2, 3)  drive pair of each axis j: secondary (1) and primary
  !>                                    (2) global node ids (periodic)
  !>   drive_loc       integer(ik)(2, 3)  their local node when owned here, else 0
  !>   drive_factor    real(rk)(3)      sigma_j / l_j of each axis (sign of the drive vector
  !>                                    over the length) [1/mm]
  type :: constraint_t
    logical :: active = .false.
    integer(ik) :: num_rows_glob = 0_ik
    integer(ik), allocatable :: glob_dof(:)
    integer(ik), allocatable :: glob_start(:)
    integer(ik), allocatable :: glob_col(:)
    real(rk), allocatable :: glob_coef(:)
    real(rk), allocatable :: glob_offset(:)
    integer(ik), allocatable :: all_cols(:)
    integer(ik) :: num_rows = 0_ik
    integer(ik), allocatable :: row_dof(:)
    integer(ik), allocatable :: row_start(:)
    real(rk), allocatable :: col_coef(:)
    integer(ik), allocatable :: col_owned(:)
    integer(ik), allocatable :: col_slot(:)
    integer(ik), allocatable :: col_local(:)
    integer(ik), allocatable :: col_glob(:)
    real(rk), allocatable :: glob_rate_coef(:, :)
    integer(ik), allocatable :: glob_prim(:)
    real(rk), allocatable :: row_offset(:)
    real(rk), allocatable :: row_offset0(:)
    real(rk), allocatable :: row_rate_coef(:, :)
    integer(ik), allocatable :: row_prim(:)
    integer(ik), allocatable :: row_of_dof(:)
    logical, allocatable :: is_constrained(:)
    type(constraint_exchange_t) :: xch
    logical :: diag_exact = .true.
    logical :: zero_mean = .false.
    logical :: mean_fixed(3) = .false.
    real(rk), allocatable :: node_weight(:)
    real(rk) :: volume0 = 0.0_rk
    real(rk) :: moment0(3) = 0.0_rk
    integer(ik) :: drive_node(2, 3) = 0_ik
    integer(ik) :: drive_loc(2, 3) = 0_ik
    real(rk) :: drive_factor(3) = 0.0_rk
  end type constraint_t

contains

  ! ==========================================================================================
  ! Building
  ! ==========================================================================================

  !> The replicated table of MPC lines (v2 read_mpcs and calc_mpcs; module header). Pure in
  !> effect (no communication); errors are fatal (every rank finds the same).
  !>
  !> Arguments:
  !>   specs       (in)   type(mpc_spec_t)(:)   MPC relations, in configuration order
  !>   nsets_glob  (in)   type(mesh_nset_t)(:)  whole-mesh node sets (global ids, file order)
  !>   cons        (inout) type(constraint_t)   the glob_* components and all_cols are set
  !> Locals:
  !>   e_dof, e_parent  integer(ik)(:)  secondary global dof and its (immediate) primary dof of
  !>                                    each entry
  !>   e_coef, e_off    real(rk)(:)     coefficient and offset of each entry
  !>   e_line           integer(ik)(:)  configuration line of each entry
  !>   order            integer(ik)(:)  entries sorted by secondary dof
  !>   s_dof            integer(ik)(:)  secondary dofs in increasing order
  !>   keep             logical(:)      sorted entry kept (not a repeated identical relation)
  !>   num_e            integer(ik)     entries
  !>   i_spec, i_set    integer         spec and node-set indices
  !>   spec_set         integer(:)      node set of each spec (0: fewer than two nodes)
  !>   i_node           integer(ik)     position in the node set
  !>   p, s             integer(ik)     primary and secondary node ids
  !>   r, pos, steps    integer(ik)     row, position of the current primary among the
  !>                                    secondaries, chain steps
  !>   cur              integer(ik)     current primary dof while resolving a chain
  !>   c, d             real(rk)        composed coefficient and offset
  !>   text             character(160)  message
  subroutine constraint_table_mpc(specs, nsets_glob, cons)

    type(mpc_spec_t), intent(in) :: specs(:)
    type(mesh_nset_t), intent(in) :: nsets_glob(:)
    type(constraint_t), intent(inout) :: cons

    integer(ik), allocatable :: e_dof(:), e_parent(:), e_line(:), order(:), s_dof(:)
    real(rk), allocatable :: e_coef(:), e_off(:)
    integer(ik) :: num_e
    integer :: i_spec, i_set
    integer, allocatable :: spec_set(:)
    logical, allocatable :: keep(:)
    integer(ik) :: i_node, p, s, r, pos, steps, cur
    real(rk) :: c, d
    character(len=160) :: text

    ! The node set of each spec (an error if it is not in the mesh; a warning, and no
    ! relation, with fewer than two nodes) and the number of entries.
    allocate (spec_set(size(specs)))
    num_e = 0_ik
    do i_spec = 1, size(specs)
      i_set = 0
      do i_node = 1_ik, int(size(nsets_glob), ik)
        if (nsets_glob(i_node)%label == specs(i_spec)%nset) i_set = int(i_node)
      end do
      if (i_set == 0) call fatal_error('constraint_table_mpc', 'line ' // &
        & int_text(specs(i_spec)%line) // ': node set ' // specs(i_spec)%nset // &
        & ' of a multi-point constraint is not in the mesh')
      spec_set(i_spec) = i_set
      if (size(nsets_glob(i_set)%nodes) < 2) then
        call log_warning('Line ' // int_text(specs(i_spec)%line) // ': node set ' // &
          & specs(i_spec)%nset // ' has fewer than two nodes: no multi-point constraint.', 2)
        spec_set(i_spec) = 0
        cycle
      end if
      num_e = num_e + int(count(nsets_glob(i_set)%nodes(2:) /= nsets_glob(i_set)%nodes(1)), ik)
    end do

    ! Entries: every node of each node set but the first is a secondary of the first (v2
    ! read_mpcs), for the spec's component.
    allocate (e_dof(num_e), e_parent(num_e), e_line(num_e), e_coef(num_e), e_off(num_e))
    num_e = 0_ik
    do i_spec = 1, size(specs)
      if (spec_set(i_spec) == 0) cycle
      associate (nodes => nsets_glob(spec_set(i_spec))%nodes)
        p = nodes(1)
        do i_node = 2_ik, int(size(nodes), ik)
          s = nodes(i_node)
          if (s == p) cycle
          num_e = num_e + 1_ik
          e_dof(num_e) = 3_ik*(s - 1_ik) + int(specs(i_spec)%comp, ik)
          e_parent(num_e) = 3_ik*(p - 1_ik) + int(specs(i_spec)%comp, ik)
          e_coef(num_e) = specs(i_spec)%coeff
          e_off(num_e) = specs(i_spec)%offset
          e_line(num_e) = specs(i_spec)%line
        end do
      end associate
    end do

    ! Entries by secondary dof. A dof that two lines make a secondary with the same relation
    ! (the same primary, coefficient and offset, e.g. an edge node of two faces whose node sets
    ! start with the same corner) keeps one entry; with different relations it is an error
    ! (D73).
    order = sort_index(e_dof)
    allocate (keep(num_e))
    keep = .true.
    do r = 2_ik, num_e
      if (e_dof(order(r)) /= e_dof(order(r - 1_ik))) cycle
      if (e_parent(order(r)) == e_parent(order(r - 1_ik)) .and. .not. &
        & (abs(e_coef(order(r)) - e_coef(order(r - 1_ik))) > 0.0_rk) .and. .not. &
        & (abs(e_off(order(r)) - e_off(order(r - 1_ik))) > 0.0_rk)) then
        keep(r) = .false.
        cycle
      end if
      write (text, '(a,i0,a,a,a,i0,a,i0,a)') 'node ', (e_dof(order(r)) - 1_ik)/3_ik + 1_ik, &
        & ' component ', axis_letter(mod(e_dof(order(r)) - 1_ik, 3_ik) + 1_ik), &
        & ' is a secondary of two different multi-point constraints (lines ', &
        & min(e_line(order(r - 1_ik)), e_line(order(r))), ' and ', &
        & max(e_line(order(r - 1_ik)), e_line(order(r))), ')'
      call fatal_error('constraint_table_mpc', trim(text))
    end do
    order = pack(order, keep)
    num_e = int(size(order), ik)
    s_dof = e_dof(order)

    ! Rows: each secondary resolved to its top primary (v2 calc_mpcs): v_s = c v_cur + d, and
    ! while cur is itself a secondary (v_cur = c2 v_next + d2), c <- c c2, d <- c d2 + d. A
    ! chain longer than the number of entries is a cycle (an error).
    cons%num_rows_glob = num_e
    allocate (cons%glob_dof(num_e), cons%glob_start(num_e + 1_ik), cons%glob_col(num_e), &
      & cons%glob_coef(num_e), cons%glob_offset(num_e))
    do r = 1_ik, num_e
      cur = e_parent(order(r))
      c = e_coef(order(r))
      d = e_off(order(r))
      steps = 0_ik
      do
        pos = find_sorted(s_dof, cur)
        if (pos == 0_ik) exit
        steps = steps + 1_ik
        if (steps > num_e) call fatal_error('constraint_table_mpc', 'the multi-point ' // &
          & 'constraints of line ' // int_text(e_line(order(r))) // ' form a cycle (a ' // &
          & 'primary node is a secondary of its own secondary)')
        d = c*e_off(order(pos)) + d
        c = c*e_coef(order(pos))
        cur = e_parent(order(pos))
      end do
      cons%glob_dof(r) = s_dof(r)
      cons%glob_start(r) = r
      cons%glob_col(r) = cur
      cons%glob_coef(r) = c
      cons%glob_offset(r) = d
    end do
    cons%glob_start(num_e + 1_ik) = num_e + 1_ik

    ! No rate offsets; the top primary is the row's one column.
    allocate (cons%glob_rate_coef(9, num_e))
    cons%glob_rate_coef = 0.0_rk
    cons%glob_prim = cons%glob_col

    ! Distinct columns, increasing.
    cons%all_cols = unique_sorted(cons%glob_col)

  end subroutine constraint_table_mpc

  !> The MPC table, its local rows and the exchange (module policy). Collective.
  !>
  !> Arguments:
  !>   cons   (out)  type(constraint_t)    the constraints (active if any row exists)
  !>   mesh   (in)   type(mesh_t)          this rank's mesh (nsets_glob, local numbering)
  !>   specs  (in)   type(mpc_spec_t)(:)   MPC relations
  !> Locals:
  !>   num_cols_glob  integer(ik)  columns of the whole table
  subroutine constraint_build_mpc(cons, mesh, specs)

    type(constraint_t), intent(out) :: cons
    type(mesh_t), intent(in) :: mesh
    type(mpc_spec_t), intent(in) :: specs(:)

    integer(ik) :: num_cols_glob

    ! The replicated table, then the local rows and the transfer.
    call constraint_table_mpc(specs, mesh%nsets_glob, cons)
    cons%active = cons%num_rows_glob > 0_ik
    call constraint_localize(cons, mesh)
    num_cols_glob = int(size(cons%all_cols), ik)
    call log_info('Multi-point constraints: ' // int_text(cons%num_rows_glob) // &
      & ' constrained dofs, ' // int_text(num_cols_glob) // ' primary dofs', 2)

  end subroutine constraint_build_mpc

  !> Local rows of the replicated table and the exchange set-up (module policy); the offsets
  !> are the constant ones (constraint_set_rates adds those of the rates). Collective.
  !>
  !> Arguments:
  !>   cons  (inout)  type(constraint_t)  table built (glob_* and all_cols); local components
  !>                                      set
  !>   mesh  (in)     type(mesh_t)        this rank's mesh
  !> Locals:
  !>   needed       logical(:)      columns referenced by the local rows
  !>   slot_of_col  integer(ik)(:)  exchange slot of each column of all_cols (0: none)
  !>   r, e, n_e    integer(ik)     global row, column entry, local entries
  !>   i_loc, k     integer(ik)     local node and component
  !>   pos          integer(ik)     position in all_cols
  !>   n_rows       integer(ik)     local rows placed
  !>   n_owned      integer(ik)     owned dofs
  subroutine constraint_localize(cons, mesh)

    type(constraint_t), intent(inout) :: cons
    type(mesh_t), intent(in) :: mesh

    logical, allocatable :: needed(:)
    integer(ik), allocatable :: slot_of_col(:)
    integer(ik) :: r, e, n_e, i_loc, k, pos, n_rows, n_owned

    n_owned = 3_ik*mesh%num_nodes_owned
    allocate (cons%is_constrained(3_ik*mesh%num_nodes_loc), &
      & cons%row_of_dof(3_ik*mesh%num_nodes_loc))
    cons%is_constrained = .false.
    cons%row_of_dof = 0_ik

    ! Local constrained dofs (owned and ghost) and the columns they need.
    allocate (needed(size(cons%all_cols)))
    needed = .false.
    n_rows = 0_ik
    n_e = 0_ik
    do r = 1_ik, cons%num_rows_glob
      i_loc = node_loc_from_glob(mesh, (cons%glob_dof(r) - 1_ik)/3_ik + 1_ik)
      if (i_loc == 0_ik) cycle
      k = mod(cons%glob_dof(r) - 1_ik, 3_ik) + 1_ik
      n_rows = n_rows + 1_ik
      cons%is_constrained(3_ik*(i_loc - 1_ik) + k) = .true.
      do e = cons%glob_start(r), cons%glob_start(r + 1_ik) - 1_ik
        n_e = n_e + 1_ik
        needed(find_sorted(cons%all_cols, cons%glob_col(e))) = .true.
      end do
    end do

    ! The transfer of the columns not owned here (collective).
    call constraint_exchange_setup(cons%xch, mesh, cons%all_cols, needed)
    allocate (slot_of_col(size(cons%all_cols)))
    slot_of_col = 0_ik
    do e = 1_ik, cons%xch%num_slots
      slot_of_col(find_sorted(cons%all_cols, cons%xch%slot_dof_glob(e))) = e
    end do

    ! Local rows, in increasing local dof, with their columns.
    cons%num_rows = n_rows
    allocate (cons%row_dof(n_rows), cons%row_start(n_rows + 1_ik), cons%row_offset(n_rows))
    allocate (cons%row_offset0(n_rows), cons%row_rate_coef(9, n_rows), cons%row_prim(n_rows))
    allocate (cons%col_coef(n_e), cons%col_owned(n_e), cons%col_slot(n_e), &
      & cons%col_local(n_e), cons%col_glob(n_e))
    cons%row_dof = pack([(r, r = 1_ik, 3_ik*mesh%num_nodes_loc)], cons%is_constrained)
    n_e = 0_ik
    do n_rows = 1_ik, cons%num_rows
      ! The global row of this local dof.
      i_loc = (cons%row_dof(n_rows) - 1_ik)/3_ik + 1_ik
      k = mod(cons%row_dof(n_rows) - 1_ik, 3_ik) + 1_ik
      r = find_sorted(cons%glob_dof, 3_ik*(mesh%node_id_glob(i_loc) - 1_ik) + k)
      cons%row_of_dof(cons%row_dof(n_rows)) = n_rows
      cons%row_start(n_rows) = n_e + 1_ik
      cons%row_offset(n_rows) = cons%glob_offset(r)
      cons%row_offset0(n_rows) = cons%glob_offset(r)
      cons%row_rate_coef(:, n_rows) = cons%glob_rate_coef(:, r)
      cons%row_prim(n_rows) = 0_ik
      do e = cons%glob_start(r), cons%glob_start(r + 1_ik) - 1_ik
        n_e = n_e + 1_ik
        if (cons%glob_col(e) == cons%glob_prim(r)) cons%row_prim(n_rows) = n_e
        cons%col_coef(n_e) = cons%glob_coef(e)
        cons%col_glob(n_e) = cons%glob_col(e)
        i_loc = node_loc_from_glob(mesh, (cons%glob_col(e) - 1_ik)/3_ik + 1_ik)
        k = mod(cons%glob_col(e) - 1_ik, 3_ik) + 1_ik
        cons%col_local(n_e) = 0_ik
        if (i_loc > 0_ik) cons%col_local(n_e) = 3_ik*(i_loc - 1_ik) + k
        cons%col_owned(n_e) = 0_ik
        if (cons%col_local(n_e) >= 1_ik .and. cons%col_local(n_e) <= n_owned) &
          & cons%col_owned(n_e) = cons%col_local(n_e)
        pos = find_sorted(cons%all_cols, cons%glob_col(e))
        cons%col_slot(n_e) = 0_ik
        if (cons%col_owned(n_e) == 0_ik) cons%col_slot(n_e) = slot_of_col(pos)
        if (cons%col_owned(n_e) == 0_ik .and. cons%col_slot(n_e) == 0_ik) call internal_error( &
          & 'constraint_build', 'a constraint column is neither owned nor a slot')
      end do
      if (cons%row_prim(n_rows) == 0_ik) call internal_error('constraint_build', &
        & 'the top primary of a constraint row is not one of its columns')
    end do
    cons%row_start(cons%num_rows + 1_ik) = n_e + 1_ik

  end subroutine constraint_localize

  ! ==========================================================================================
  ! Operations
  ! ==========================================================================================

  !> v = T v_hat (+ g) on the local constrained dofs (module policy). Collective.
  !>
  !> Arguments:
  !>   cons          (in)     type(constraint_t)  the constraints
  !>   x             (inout)  real(rk)(:)  (3 num_nodes_loc) in: valid at the owned dofs (the
  !>                                       independent values; constrained entries ignored);
  !>                                       out: every local constrained dof set
  !>   with_offsets  (in)     logical      add the offsets g (a velocity) or not (a correction)
  !> Locals:
  !>   buf    real(rk)(:)  values of the slots
  !>   r, e   integer(ik)  local row and column entry
  !>   value  real(rk)     the row's value
  subroutine constraint_expand(cons, x, with_offsets)

    type(constraint_t), intent(in) :: cons
    real(rk), intent(inout) :: x(:)
    logical, intent(in) :: with_offsets

    real(rk), allocatable :: buf(:)
    integer(ik) :: r, e
    real(rk) :: value

    ! The columns owned elsewhere.
    allocate (buf(cons%xch%num_slots))
    call constraint_forward(cons%xch, x, buf)

    ! Each local constrained dof from its columns (never constrained themselves).
    do r = 1_ik, cons%num_rows
      value = 0.0_rk
      if (with_offsets) value = cons%row_offset(r)
      do e = cons%row_start(r), cons%row_start(r + 1_ik) - 1_ik
        if (cons%col_owned(e) > 0_ik) then
          value = value + cons%col_coef(e)*x(cons%col_owned(e))
        else
          value = value + cons%col_coef(e)*buf(cons%col_slot(e))
        end if
      end do
      x(cons%row_dof(r)) = value
    end do

  end subroutine constraint_expand

  !> y = T^T y on the independent dofs, zero at the constrained dofs (module policy).
  !> Collective.
  !>
  !> Arguments:
  !>   cons     (in)     type(constraint_t)  the constraints
  !>   n_owned  (in)     integer(ik)         owned dofs (3 num_nodes_owned)
  !>   y        (inout)  real(rk)(:)  (3 num_nodes_loc) in: complete at the owned dofs; out:
  !>                                  complete at the owned independent dofs, zero at every
  !>                                  local constrained dof (other ghost entries unchanged)
  !> Locals:
  !>   buf   real(rk)(:)  contributions to the slots
  !>   r, e  integer(ik)  local row and column entry
  !>   ys    real(rk)     value at the constrained dof
  subroutine constraint_reduce(cons, n_owned, y)

    type(constraint_t), intent(in) :: cons
    integer(ik), intent(in) :: n_owned
    real(rk), intent(inout) :: y(:)

    real(rk), allocatable :: buf(:)
    integer(ik) :: r, e
    real(rk) :: ys

    ! Contributions of the owned constrained dofs to their columns.
    allocate (buf(cons%xch%num_slots))
    buf = 0.0_rk
    do r = 1_ik, cons%num_rows
      if (cons%row_dof(r) > n_owned) cycle
      ys = y(cons%row_dof(r))
      do e = cons%row_start(r), cons%row_start(r + 1_ik) - 1_ik
        if (cons%col_owned(e) > 0_ik) then
          y(cons%col_owned(e)) = y(cons%col_owned(e)) + cons%col_coef(e)*ys
        else
          buf(cons%col_slot(e)) = buf(cons%col_slot(e)) + cons%col_coef(e)*ys
        end if
      end do
    end do

    ! Contributions to columns owned elsewhere; zero at the constrained dofs.
    call constraint_reverse_add(cons%xch, buf, y)
    y(cons%row_dof) = 0.0_rk

  end subroutine constraint_reduce

  !> y = T^T K T x by elements (R-P4 Algorithm 1: T, the element product, T^T), zero at the
  !> constrained and the prescribed dofs. Collective.
  !>
  !> Arguments:
  !>   cons     (in)   type(constraint_t)  the constraints
  !>   mesh     (in)   type(mesh_t)        rank-local mesh
  !>   k_elt    (in)   real(rk)(:, :, :)   element matrices
  !>   bc_mask  (in)   logical(:)          (3 num_nodes_loc) prescribed dofs
  !>   x        (in)   real(rk)(:)         (3 num_nodes_loc) reduced vector, valid on the owned
  !>                                       and the unconstrained ghost dofs
  !>   y        (out)  real(rk)(:)         (3 num_nodes_loc) result (complete at the owned dofs)
  !> Locals:
  !>   x_full   real(rk)(:)  T x
  subroutine constraint_matvec(cons, mesh, k_elt, bc_mask, x, y)

    type(constraint_t), intent(in) :: cons
    type(mesh_t), intent(in) :: mesh
    real(rk), intent(in) :: k_elt(:, :, :)
    logical, intent(in) :: bc_mask(:)
    real(rk), intent(in) :: x(:)
    real(rk), intent(out) :: y(:)

    real(rk), allocatable :: x_full(:)

    ! T x, the element product (zero at the prescribed dofs), T^T, zero at the prescribed dofs.
    x_full = x
    call constraint_expand(cons, x_full, .false.)
    call ebe_matvec(mesh, k_elt, bc_mask, x_full, y)
    call constraint_reduce(cons, 3_ik*mesh%num_nodes_owned, y)
    where (bc_mask) y = 0.0_rk

  end subroutine constraint_matvec

  !> Exact diagonal of T^T K T: for each independent dof q, the sum over the elements of
  !> sum_{a,b in e} T_aq K_ab T_bq (v2-fixed assemble_diagonals("MPC") after C67, addendum §1.1,
  !> §3.1), complete on the owned and the ghost dofs, zero at the constrained dofs. Collective.
  !> (The lumped periodic diagonal: constraint_diagonal_lumped; the caller chooses by diag_exact.)
  !>
  !> Arguments:
  !>   cons   (in)   type(constraint_t)  the constraints
  !>   mesh   (in)   type(mesh_t)        rank-local mesh
  !>   k_elt  (in)   real(rk)(:, :, :)   element matrices
  !>   diag   (out)  real(rk)(:)         (3 num_nodes_loc) the diagonal
  !> Locals:
  !>   buf         real(rk)(:)     contributions to the slots
  !>   n           integer(ik)     dofs per element
  !>   ent_pos     integer(ik)(:)  element dof of each entry of the element's expansion
  !>   ent_key     integer(ik)(:)  independent global dof of each entry
  !>   ent_coef    real(rk)(:)     coefficient of each entry
  !>   ent_local   integer(ik)(:)  local dof of the entry's independent dof (0: not local)
  !>   ent_slot    integer(ik)(:)  slot of the entry's independent dof (when not local)
  !>   n_ent       integer(ik)     entries of the element
  !>   i_elt, i, a, b  integer(ik) element, element dof, entries
  !>   dof, r, e   integer(ik)     local dof, local row, column entry
  !>   contrib     real(rk)        contribution to an independent dof
  !>   has_cons    logical         the element has a constrained dof
  !>   max_cols    integer(ik)     most columns of a local row (entries per element dof)
  subroutine constraint_diagonal(cons, mesh, k_elt, diag)

    type(constraint_t), intent(in) :: cons
    type(mesh_t), intent(in) :: mesh
    real(rk), intent(in) :: k_elt(:, :, :)
    real(rk), intent(out) :: diag(:)

    real(rk), allocatable :: buf(:)
    integer(ik) :: n
    integer(ik), allocatable :: ent_pos(:), ent_key(:), ent_local(:), ent_slot(:)
    real(rk), allocatable :: ent_coef(:)
    integer(ik) :: n_ent, i_elt, i, a, b, dof, r, e, max_cols
    real(rk) :: contrib
    logical :: has_cons

    n = int(size(k_elt, 1), ik)
    max_cols = 1_ik
    do r = 1_ik, cons%num_rows
      max_cols = max(max_cols, cons%row_start(r + 1_ik) - cons%row_start(r))
    end do
    allocate (ent_pos(n*max_cols), ent_key(n*max_cols), ent_coef(n*max_cols), &
      & ent_local(n*max_cols), ent_slot(n*max_cols))
    allocate (buf(cons%xch%num_slots))
    buf = 0.0_rk
    diag = 0.0_rk
    do i_elt = 1_ik, mesh%num_elts
      ! Elements without constrained dofs: their diagonal (as ebe_diagonal).
      has_cons = .false.
      do i = 1_ik, n
        if (cons%is_constrained(mesh%elt_dofs_loc(i, i_elt))) has_cons = .true.
      end do
      if (.not. has_cons) then
        do i = 1_ik, n
          dof = mesh%elt_dofs_loc(i, i_elt)
          diag(dof) = diag(dof) + k_elt(i, i, i_elt)
        end do
        cycle
      end if

      ! Expansion of the element dofs into independent dofs: an unconstrained dof is itself
      ! (coefficient 1), a constrained dof its row's columns.
      n_ent = 0_ik
      do i = 1_ik, n
        dof = mesh%elt_dofs_loc(i, i_elt)
        r = cons%row_of_dof(dof)
        if (r == 0_ik) then
          n_ent = n_ent + 1_ik
          ent_pos(n_ent) = i
          ent_key(n_ent) = 3_ik*(mesh%node_id_glob((dof - 1_ik)/3_ik + 1_ik) - 1_ik) + &
            & mod(dof - 1_ik, 3_ik) + 1_ik
          ent_coef(n_ent) = 1.0_rk
          ent_local(n_ent) = dof
          ent_slot(n_ent) = 0_ik
        else
          do e = cons%row_start(r), cons%row_start(r + 1_ik) - 1_ik
            n_ent = n_ent + 1_ik
            ent_pos(n_ent) = i
            ent_key(n_ent) = cons%col_glob(e)
            ent_coef(n_ent) = cons%col_coef(e)
            ent_local(n_ent) = cons%col_local(e)
            ent_slot(n_ent) = cons%col_slot(e)
          end do
        end if
      end do

      ! For each entry a, its independent dof q gets sum over the entries b with the same q of
      ! T_aq K_ab T_bq (each entry a once: the full double sum over a and b).
      do a = 1_ik, n_ent
        contrib = 0.0_rk
        do b = 1_ik, n_ent
          if (ent_key(b) /= ent_key(a)) cycle
          contrib = contrib + ent_coef(a)*k_elt(ent_pos(a), ent_pos(b), i_elt)*ent_coef(b)
        end do
        if (ent_local(a) > 0_ik) then
          diag(ent_local(a)) = diag(ent_local(a)) + contrib
        else
          buf(ent_slot(a)) = buf(ent_slot(a)) + contrib
        end if
      end do
    end do

    ! Ghost and slot contributions to their owners; zero at the constrained dofs; the ghosts.
    call halo_reverse_add(mesh, diag, 3_ik)
    call constraint_reverse_add(cons%xch, buf, diag)
    where (cons%is_constrained) diag = 0.0_rk
    call halo_forward(mesh, diag, 3_ik)

  end subroutine constraint_diagonal

  !> Velocity conditions carried through the MPC groups (v2-fixed propagate_bcs_mpcs, C66;
  !> addendum §3.2): each condition on a node of a group gives its primary (v_given - d)/c, the
  !> primary's own condition v_given; the values of a group are averaged; every dof of the group
  !> is then prescribed at c v_p + d. A given condition that differs from it by more than
  !> 1e-10 of the largest condition is an error. Collective.
  !>
  !> Arguments:
  !>   cons     (in)     type(constraint_t)  the constraints (MPCs: one column per row)
  !>   mesh     (in)     type(mesh_t)        rank-local mesh
  !>   bc_mask  (inout)  logical(:)          (3 num_nodes_loc) prescribed dofs (owned and ghost)
  !>   bc_vel   (inout)  real(rk)(:)         (3 num_nodes_loc) prescribed values [mm/s]
  !> Locals:
  !>   n_owned     integer(ik)  owned dofs
  !>   given_mask  logical(:)   the conditions as given
  !>   given_vel   real(rk)(:)  their values [mm/s]
  !>   sum_v, cnt  real(rk)(:)  per independent dof: sum and number of the candidate values
  !>   b_sum, b_cnt  real(rk)(:)  the same for the slots
  !>   val, flag   real(rk)(:)  per independent dof: the group value and whether it is set
  !>   b_val, b_flag  real(rk)(:)  the same for the slots
  !>   r, e, s, i  integer(ik)  local row, column entry, constrained dof, dof
  !>   cand        real(rk)     candidate value of a primary [mm/s]
  !>   pv, pf      real(rk)     the primary's value and flag
  !>   scale       real(rk)     largest |given condition| over all ranks [mm/s]
  !>   bad, any_bad  logical    a conflict on this rank, on any rank
  !>   text        character(200) message
  subroutine constraint_propagate_bcs(cons, mesh, bc_mask, bc_vel)

    type(constraint_t), intent(in) :: cons
    type(mesh_t), intent(in) :: mesh
    logical, intent(inout) :: bc_mask(:)
    real(rk), intent(inout) :: bc_vel(:)

    integer(ik) :: n_owned
    logical, allocatable :: given_mask(:)
    real(rk), allocatable :: given_vel(:), sum_v(:), cnt(:), b_sum(:), b_cnt(:)
    real(rk), allocatable :: val(:), flag(:), b_val(:), b_flag(:)
    integer(ik) :: r, e, s, i
    real(rk) :: cand, pv, pf, scale
    logical :: bad, any_bad
    character(len=200) :: text

    n_owned = 3_ik*mesh%num_nodes_owned
    given_mask = bc_mask
    given_vel = bc_vel
    allocate (sum_v(size(bc_vel)), cnt(size(bc_vel)), val(size(bc_vel)), flag(size(bc_vel)))
    allocate (b_sum(cons%xch%num_slots), b_cnt(cons%xch%num_slots))
    allocate (b_val(cons%xch%num_slots), b_flag(cons%xch%num_slots))
    sum_v = 0.0_rk
    cnt = 0.0_rk
    b_sum = 0.0_rk
    b_cnt = 0.0_rk

    ! Candidates: an owned unconstrained dof with a condition is its own; an owned constrained
    ! dof with a condition gives its primary (v - d)/c (c is nonzero: config_mod).
    do i = 1_ik, n_owned
      if (.not. given_mask(i) .or. cons%is_constrained(i)) cycle
      sum_v(i) = sum_v(i) + given_vel(i)
      cnt(i) = cnt(i) + 1.0_rk
    end do
    do r = 1_ik, cons%num_rows
      s = cons%row_dof(r)
      if (s > n_owned .or. .not. given_mask(s)) cycle
      e = cons%row_start(r)
      cand = (given_vel(s) - cons%row_offset(r))/cons%col_coef(e)
      if (cons%col_owned(e) > 0_ik) then
        sum_v(cons%col_owned(e)) = sum_v(cons%col_owned(e)) + cand
        cnt(cons%col_owned(e)) = cnt(cons%col_owned(e)) + 1.0_rk
      else
        b_sum(cons%col_slot(e)) = b_sum(cons%col_slot(e)) + cand
        b_cnt(cons%col_slot(e)) = b_cnt(cons%col_slot(e)) + 1.0_rk
      end if
    end do
    call constraint_reverse_add(cons%xch, b_sum, sum_v)
    call constraint_reverse_add(cons%xch, b_cnt, cnt)

    ! The averaged value of each owned independent dof with candidates; the ghosts and slots.
    val = 0.0_rk
    flag = 0.0_rk
    do i = 1_ik, n_owned
      if (cons%is_constrained(i) .or. cnt(i) < 0.5_rk) cycle
      val(i) = sum_v(i)/cnt(i)
      flag(i) = 1.0_rk
    end do
    call halo_forward(mesh, val, 3_ik)
    call halo_forward(mesh, flag, 3_ik)
    call constraint_forward(cons%xch, val, b_val)
    call constraint_forward(cons%xch, flag, b_flag)

    ! New conditions: the independent dofs from their values, the constrained dofs from their
    ! primary (c v_p + d) when it is set.
    do i = 1_ik, int(size(bc_mask), ik)
      if (cons%is_constrained(i)) cycle
      bc_mask(i) = flag(i) > 0.5_rk
      if (bc_mask(i)) bc_vel(i) = val(i)
    end do
    do r = 1_ik, cons%num_rows
      s = cons%row_dof(r)
      e = cons%row_start(r)
      if (cons%col_owned(e) > 0_ik) then
        pv = val(cons%col_owned(e))
        pf = flag(cons%col_owned(e))
      else
        pv = b_val(cons%col_slot(e))
        pf = b_flag(cons%col_slot(e))
      end if
      bc_mask(s) = pf > 0.5_rk
      if (bc_mask(s)) bc_vel(s) = cons%col_coef(e)*pv + cons%row_offset(r)
    end do

    ! Conflicts: a given condition that the group does not reproduce (C66).
    scale = 0.0_rk
    do i = 1_ik, n_owned
      if (given_mask(i)) scale = max(scale, abs(given_vel(i)))
    end do
    call par_max(scale, pv)
    scale = pv
    bad = .false.
    do i = 1_ik, n_owned
      if (.not. given_mask(i)) cycle
      if (abs(bc_vel(i) - given_vel(i)) > conflict_rel_tol*scale) then
        if (.not. bad) then
          write (text, '(a,i0,a,a,a,es12.5,a,es12.5,a)') 'Node ', &
            & mesh%node_id_glob((i - 1_ik)/3_ik + 1_ik), ' component ', &
            & axis_letter(mod(i - 1_ik, 3_ik) + 1_ik), ': condition ', given_vel(i), &
            & ', from its multi-point constraint group ', bc_vel(i), '.'
          call log_warning(trim(text), 1, all_ranks=.true.)
        end if
        bad = .true.
      end if
    end do
    ! (Every rank's warning written by rank 0 before the collective error, step 9.7.)
    call par_log_flush()
    call par_any(bad, any_bad)
    if (any_bad) call fatal_error('constraint_propagate_bcs', 'conflicting velocity ' // &
      & 'conditions in a multi-point constraint group')

  end subroutine constraint_propagate_bcs

  !> A velocity projected onto the constraints: every local constrained dof set to
  !> sum T_sq v_q + g_s from the independent values (v2's uniaxial driver after its initial
  !> estimate, driver_uniaxial_control_mod1.f90 l. 94-103). Collective.
  !>
  !> Arguments:
  !>   cons  (in)     type(constraint_t)  the constraints
  !>   vel   (inout)  real(rk)(:)         (3 num_nodes_loc) velocities, valid at the owned and
  !>                                      the unconstrained ghost dofs [mm/s]
  subroutine constraint_project(cons, vel)

    type(constraint_t), intent(in) :: cons
    real(rk), intent(inout) :: vel(:)

    call constraint_expand(cons, vel, .true.)

  end subroutine constraint_project

  !> The lumped periodic Jacobi diagonal (periodic_preconditioner lumped_diagonal, the default;
  !> v2-fixed conjugate_gradient_mod1.f90::assemble_diagonals
  !> "PBC" with pre_cond 1; R-P4 Algorithm 1 l. 5, M_hat = (T+)^T M, addendum §1.2, §3.1): every
  !> element diagonal K_ii is added to the dof itself when it is independent, otherwise to the
  !> top primary of its row (row_prim); the drive columns are left out. Complete on the owned and
  !> the ghost dofs, zero at the constrained dofs. Collective.
  !>
  !> Arguments:
  !>   cons   (in)   type(constraint_t)  the constraints
  !>   mesh   (in)   type(mesh_t)        rank-local mesh
  !>   k_elt  (in)   real(rk)(:, :, :)   element matrices
  !>   diag   (out)  real(rk)(:)         (3 num_nodes_loc) the diagonal
  !> Locals:
  !>   buf        real(rk)(:)  contributions to the slots
  !>   i_elt, i   integer(ik)  element, element dof
  !>   dof, r, e  integer(ik)  local dof, local row, its top-primary column entry
  subroutine constraint_diagonal_lumped(cons, mesh, k_elt, diag)

    type(constraint_t), intent(in) :: cons
    type(mesh_t), intent(in) :: mesh
    real(rk), intent(in) :: k_elt(:, :, :)
    real(rk), intent(out) :: diag(:)

    real(rk), allocatable :: buf(:)
    integer(ik) :: i_elt, i, dof, r, e

    allocate (buf(cons%xch%num_slots))
    buf = 0.0_rk
    diag = 0.0_rk
    do i_elt = 1_ik, mesh%num_elts
      do i = 1_ik, int(size(k_elt, 1), ik)
        dof = mesh%elt_dofs_loc(i, i_elt)
        r = cons%row_of_dof(dof)
        if (r == 0_ik) then
          ! An independent dof: its own diagonal.
          diag(dof) = diag(dof) + k_elt(i, i, i_elt)
        else
          ! A constrained dof: onto its top primary (local, or through the exchange).
          e = cons%row_prim(r)
          if (cons%col_local(e) > 0_ik) then
            diag(cons%col_local(e)) = diag(cons%col_local(e)) + k_elt(i, i, i_elt)
          else
            buf(cons%col_slot(e)) = buf(cons%col_slot(e)) + k_elt(i, i, i_elt)
          end if
        end if
      end do
    end do

    ! Ghost and slot contributions to their owners; zero at the constrained dofs; the ghosts.
    call halo_reverse_add(mesh, diag, 3_ik)
    call constraint_reverse_add(cons%xch, buf, diag)
    where (cons%is_constrained) diag = 0.0_rk
    call halo_forward(mesh, diag, 3_ik)

  end subroutine constraint_diagonal_lumped

  !> Offsets of the imposed macroscopic rates of a step (periodic conditions, step 8.3b):
  !> g_s = offset0_s + sum_c a_sc rate_c on every local row.
  !>
  !> Arguments:
  !>   cons   (inout)  type(constraint_t)  the constraints (row_offset set)
  !>   rates  (in)     real(rk)(9)         imposed rates of Fbar-dot, index 3 (i - 1) + j (0
  !>                                       where free) [1/s]
  !> Locals:
  !>   r      integer(ik)  local row
  subroutine constraint_set_rates(cons, rates)

    type(constraint_t), intent(inout) :: cons
    real(rk), intent(in) :: rates(9)

    integer(ik) :: r

    do r = 1_ik, cons%num_rows
      cons%row_offset(r) = cons%row_offset0(r) + dot_product(cons%row_rate_coef(:, r), rates)
    end do

  end subroutine constraint_set_rates

  !> Velocity conditions with periodic conditions (Q-P6, addendum §3.3): applied as given (v2,
  !> no propagation); a condition on a constrained dof cannot hold through v = T v_hat + g and
  !> is an error; the components with a condition are fixed by it (mean_fixed: not shifted by
  !> constraint_zero_mean). Collective.
  !>
  !> Arguments:
  !>   cons     (inout)  type(constraint_t)  the constraints (mean_fixed set)
  !>   mesh     (in)     type(mesh_t)        rank-local mesh
  !>   bc_mask  (in)     logical(:)          (3 num_nodes_loc) prescribed dofs
  !> Locals:
  !>   n_owned  integer(ik)   owned dofs
  !>   i, k     integer(ik)   dof, component
  !>   bad, any_bad, part  logical  a condition on a constrained dof here, anywhere; a
  !>                                component with a condition here
  !>   text     character(200)  message
  subroutine constraint_check_bcs(cons, mesh, bc_mask)

    type(constraint_t), intent(inout) :: cons
    type(mesh_t), intent(in) :: mesh
    logical, intent(in) :: bc_mask(:)

    integer(ik) :: n_owned, i, k
    logical :: bad, any_bad, part
    character(len=200) :: text

    n_owned = 3_ik*mesh%num_nodes_owned

    ! No condition on a constrained (owned) dof.
    bad = .false.
    do i = 1_ik, n_owned
      if (.not. (bc_mask(i) .and. cons%is_constrained(i))) cycle
      if (.not. bad) then
        write (text, '(a,i0,a,a,a)') 'Node ', mesh%node_id_glob((i - 1_ik)/3_ik + 1_ik), &
          & ' component ', axis_letter(mod(i - 1_ik, 3_ik) + 1_ik), ': a velocity ' // &
          & 'condition on a dof that the periodicity relations constrain.'
        call log_warning(trim(text), 1, all_ranks=.true.)
      end if
      bad = .true.
    end do
    ! (Every rank's warning written by rank 0 before the collective error, step 9.7.)
    call par_log_flush()
    call par_any(bad, any_bad)
    if (any_bad) call fatal_error('constraint_check_bcs', 'a velocity condition is on a ' // &
      & 'periodic secondary dof (or on the drive dof of an imposed component), which the ' // &
      & 'periodicity relations constrain: give it on an independent dof (an interior or ' // &
      & 'primary node)')

    ! The components fixed by conditions.
    do k = 1_ik, 3_ik
      part = any(bc_mask(k:n_owned:3))
      call par_any(part, cons%mean_fixed(k))
    end do

  end subroutine constraint_check_bcs

  !> The rate of the macroscopic deformation gradient from the drive pairs (periodic conditions;
  !> R-P4 Eq. 72 in rate form, addendum §1.2 Eq. (3)): Fbar-dot_kj = sigma_j (v_k(Sd_j) -
  !> v_k(Pd_j)) / l_j; applied to the current coordinates, the same formula gives Fbar itself
  !> (R-P4 Eq. 72: the pair's current separation over its reference one, sigma_j l_j e_j).
  !> Collective (the owners of the drive nodes contribute their values).
  !>
  !> Arguments:
  !>   cons  (in)  type(constraint_t)  the constraints (drive pairs set)
  !>   vel   (in)  real(rk)(:)         (3 num_nodes_loc) velocities [mm/s] (or current
  !>                                   coordinates [mm]), valid at the owned dofs
  !>   fdot  (out, result) real(rk)(3, 3)  Fbar-dot_kj [1/s] (or Fbar_kj [-])
  !> Locals:
  !>   part, whole  real(rk)(18)  velocities of the drive nodes, (component, end, axis): this
  !>                              rank's owned ones, and over the ranks
  !>   j, end_pair  integer       axis, end of the pair (1 secondary, 2 primary)
  !>   i_loc        integer(ik)   local node
  function constraint_macro_rate(cons, vel) result(fdot)

    type(constraint_t), intent(in) :: cons
    real(rk), intent(in) :: vel(:)
    real(rk) :: fdot(3, 3)

    real(rk) :: part(3, 2, 3), whole(3, 2, 3)
    integer :: j, end_pair
    integer(ik) :: i_loc

    ! The drive nodes' velocities from their owners.
    part = 0.0_rk
    do j = 1, 3
      do end_pair = 1, 2
        i_loc = cons%drive_loc(end_pair, j)
        if (i_loc > 0_ik) part(:, end_pair, j) = vel(3_ik*(i_loc - 1_ik) + 1_ik: &
          & 3_ik*i_loc)
      end do
    end do
    call par_sum_drive(part, whole)

    ! Column j: the pair's velocity difference over sigma_j l_j.
    do j = 1, 3
      fdot(:, j) = cons%drive_factor(j)*(whole(:, 1, j) - whole(:, 2, j))
    end do

  end function constraint_macro_rate

  !> The velocity shifted by a uniform translation so that its fluctuation has zero weighted
  !> mean (periodic conditions, Q-P3, addendum §3.3; R-P4 below Eq. 70 in rate form): for each
  !> component k not fixed by a velocity condition, t_k = (sum_a w_a v_a,k - sum_j Fbar-dot_kj
  !> W_j) / V_0 and v_a,k <- v_a,k - t_k at every local node. A uniform translation changes
  !> neither the relations (differences of velocities) nor Fbar-dot. Collective.
  !>
  !> Arguments:
  !>   cons   (in)     type(constraint_t)  the constraints (weights and drive pairs set)
  !>   mesh   (in)     type(mesh_t)        rank-local mesh
  !>   vel    (inout)  real(rk)(:)         (3 num_nodes_loc) velocities, valid at the owned and
  !>                                       ghost dofs [mm/s]
  !> Locals:
  !>   fdot    real(rk)(3, 3)  Fbar-dot of vel [1/s]
  !>   part, whole  real(rk)(3)  sum_a w_a v_a,k over the owned nodes, over the ranks [mm^4/s]
  !>   shift   real(rk)(3)     the translation t [mm/s]
  !>   i_loc   integer(ik)     local node
  !>   k       integer         component
  subroutine constraint_zero_mean(cons, mesh, vel)

    type(constraint_t), intent(in) :: cons
    type(mesh_t), intent(in) :: mesh
    real(rk), intent(inout) :: vel(:)

    real(rk) :: fdot(3, 3), part(3), whole(3), shift(3)
    integer(ik) :: i_loc
    integer :: k

    if (.not. cons%zero_mean) return

    ! Weighted sums of the owned velocities, and Fbar-dot.
    part = 0.0_rk
    do i_loc = 1_ik, mesh%num_nodes_owned
      part = part + cons%node_weight(i_loc)*vel(3_ik*(i_loc - 1_ik) + 1_ik:3_ik*i_loc)
    end do
    call par_sum(part, whole)
    fdot = constraint_macro_rate(cons, vel)

    ! The translation of each free component, removed at every local node.
    do k = 1, 3
      shift(k) = 0.0_rk
      if (.not. cons%mean_fixed(k)) shift(k) = (whole(k) - dot_product(fdot(k, :), &
        & cons%moment0))/cons%volume0
    end do
    do i_loc = 1_ik, mesh%num_nodes_loc
      vel(3_ik*(i_loc - 1_ik) + 1_ik:3_ik*i_loc) = vel(3_ik*(i_loc - 1_ik) + 1_ik:3_ik*i_loc) - &
        & shift
    end do

  end subroutine constraint_zero_mean

  !> Sum over the ranks of the drive nodes' values (par_sum of the flattened array).
  !>
  !> Arguments:
  !>   part   (in)   real(rk)(3, 2, 3)  this rank's values
  !>   whole  (out)  real(rk)(3, 2, 3)  the sums
  !> Locals:
  !>   flat   real(rk)(18)  the sums, flattened
  subroutine par_sum_drive(part, whole)

    real(rk), intent(in) :: part(3, 2, 3)
    real(rk), intent(out) :: whole(3, 2, 3)

    real(rk) :: flat(18)

    call par_sum(reshape(part, [18]), flat)
    whole = reshape(flat, [3, 2, 3])

  end subroutine par_sum_drive

  ! ==========================================================================================
  ! Helpers
  ! ==========================================================================================

  !> Indices that sort an integer array in increasing order (stable merge sort).
  !>
  !> Arguments:
  !>   keys   (in)  integer(ik)(:)  values
  !>   order  (out, result) integer(ik)(size(keys))  keys(order) is increasing
  !> Locals:
  !>   tmp     integer(ik)(:)  merge buffer
  !>   width   integer(ik)     length of the sorted runs
  !>   lo, mid, hi  integer(ik)  bounds of the two runs merged
  !>   i, j, k integer(ik)     positions in the runs and the buffer
  !>   n       integer(ik)     number of keys
  pure function sort_index(keys) result(order)

    integer(ik), intent(in) :: keys(:)
    integer(ik) :: order(size(keys))

    integer(ik), allocatable :: tmp(:)
    integer(ik) :: width, lo, mid, hi, i, j, k, n

    n = int(size(keys), ik)
    order = [(i, i = 1_ik, n)]
    allocate (tmp(n))
    width = 1_ik
    do while (width < n)
      lo = 1_ik
      do while (lo <= n)
        ! Merge order(lo:mid) and order(mid+1:hi).
        mid = min(lo + width - 1_ik, n)
        hi = min(lo + 2_ik*width - 1_ik, n)
        i = lo
        j = mid + 1_ik
        k = lo
        do while (i <= mid .and. j <= hi)
          if (keys(order(j)) < keys(order(i))) then
            tmp(k) = order(j)
            j = j + 1_ik
          else
            tmp(k) = order(i)
            i = i + 1_ik
          end if
          k = k + 1_ik
        end do
        do while (i <= mid)
          tmp(k) = order(i)
          i = i + 1_ik
          k = k + 1_ik
        end do
        do while (j <= hi)
          tmp(k) = order(j)
          j = j + 1_ik
          k = k + 1_ik
        end do
        order(lo:hi) = tmp(lo:hi)
        lo = lo + 2_ik*width
      end do
      width = 2_ik*width
    end do

  end function sort_index

  !> The distinct values of an integer array, increasing.
  !>
  !> Arguments:
  !>   values  (in)  integer(ik)(:)  values (any order, repetitions allowed)
  !>   uniq    (out, result) integer(ik)(:), allocatable  distinct values, increasing
  !> Locals:
  !>   order   integer(ik)(:)  values(order) is increasing
  !>   i, n    integer(ik)     position, distinct values found
  !>
  !> Explicit loops (gfortran 15 -O3 warns of an unset descriptor with array expressions here,
  !> ROADMAP section G).
  pure function unique_sorted(values) result(uniq)

    integer(ik), intent(in) :: values(:)
    integer(ik), allocatable :: uniq(:)

    integer(ik), allocatable :: order(:)
    integer(ik) :: i, n

    allocate (order(size(values)))
    order = sort_index(values)
    n = 0_ik
    do i = 1_ik, int(size(values), ik)
      if (i == 1_ik) then
        n = n + 1_ik
      else if (values(order(i)) /= values(order(i - 1_ik))) then
        n = n + 1_ik
      end if
    end do
    allocate (uniq(n))
    n = 0_ik
    do i = 1_ik, int(size(values), ik)
      if (i == 1_ik) then
        n = n + 1_ik
        uniq(n) = values(order(i))
      else if (values(order(i)) /= values(order(i - 1_ik))) then
        n = n + 1_ik
        uniq(n) = values(order(i))
      end if
    end do

  end function unique_sorted

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

end module constraint_mod
