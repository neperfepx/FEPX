! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> MPI tests of periodic boundary conditions (periodic_mod, constraint_mod; ROADMAP step 8.3b;
!> dev/theory/addendum_periodic_mpc.md §3.3, §5 V8 and the unit and property tests), on 1, 2
!> and 4 ranks.
!>
!> Source : FEPX (ROADMAP step 8.3b).
!>
!> Usage: mpiexec -n N test_periodic <verification mesh dir> (tests/verification/cases/meshes).
!>
!> Meshes: periodic_box_2x1x1.msh (a Neper periodic tessellation of a 2 x 1 x 1 box: non-planar
!> "raspberry" boundaries, non-cubic) and periodic_cube.msh (a periodic one-cell cube: planar
!> faces; its three drive pairs share one primary node), distributed with partitioning "even" so
!> that relations cross the partitions; material: isotropic [s] = [C] (E = 200 GPa, nu = 0.3).
!> Cases: small strain with 11, 22, 12 imposed; deformation_gradient_rate with 11, 12, 13, 23
!> imposed (C70's configuration); deformation_gradient_rate with only 11 imposed (free rotation
!> rates: the table is still exact). Checks, for each mesh and case:
!>   - the table: one row per secondary dof but one freed drive dof per free parameter of
!>     Fbar-dot; no column is a constrained dof; each row's top primary is one of its columns;
!>   - the relations (2): v = T v_hat + g for a reduced vector v_hat (values from the global ids)
!>     and random imposed rates gives, for every relation (s, p, P) with both nodes local,
!>     v(s) - v(p) = Fbar-dot (P o l) to round-off (1e-12 of the jumps), with Fbar-dot from the
!>     drive pairs (constraint_macro_rate) equal to the imposed rates on the imposed components
!>     and symmetric with small strain (C72);
!>   - T and T^T adjoint over the owned dofs (relative 1e-12);
!>   - the reduced operator equal to the serial T^T K T x, the exact diagonal to (T e_q)^T K
!>     (T e_q) and the lumped diagonal (v2's pre_cond 1) to the serial sum of the element
!>     diagonals onto the top primaries
!>     (relative 1e-12);
!>   - the consistent nodal weights: sum w = V_0 and sum w X = int X dOmega_0 of a serial
!>     quadrature (relative 1e-12);
!>   - the zero-mean shift (constraint_zero_mean): sum w v = Fbar-dot W afterwards (1e-12);
!>   - V8: Fbar of the coordinates x = X + v (one unit of time) from the drive pairs equals the
!>     volume average of F (periodic_fbar_average, R-P4 Eqs. 68, 72) to 1e-10, and the volume
!>     average of a synthetic stress field (periodic_stress_average) equals a serial quadrature
!>     (relative 1e-12).
!>
!> Locals:
!>   ework  type(element_work_t)  element work arrays (element_work_alloc; D100)
!>   young, nu  real(rk)              Young's modulus [MPa], Poisson's ratio
!>   meshes     character(512)        verification mesh directory
!>   mesh_names character(32)(2)      the two meshes
!>   i_mesh, i_case  integer          mesh and case
!>   unit_log   integer               scratch unit receiving the log
program test_periodic

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk, ik
  use constants_mod, only: sqrt2
  use log_mod, only: log_set_unit
  use parallel_mod, only: par_init, par_finalize, par_sum, par_max
  use mesh_file_mod, only: mesh_file_t, mesh_file_load
  use mesh_mod, only: mesh_t, find_sorted, node_loc_from_glob
  use mesh_blocks_mod, only: mesh_load_blocks
  use element_mod, only: element_t, element_init, element_geometry_qpt
  use element_matrix_mod, only: element_stiffness
  use halo_mod, only: halo_forward, dist_dot
  use constraint_mod, only: constraint_t, constraint_expand, constraint_reduce, &
    & constraint_matvec, constraint_diagonal, constraint_diagonal_lumped, constraint_set_rates, &
    & constraint_macro_rate, constraint_zero_mean
  use periodic_mod, only: periodic_build, periodic_drive_pairs, periodic_stress_average, &
    & periodic_fbar_average
  use assert_mod, only: assert_true, assert_equal, assert_close, test_summary
  use element_matrix_mod, only: element_work_t, element_work_alloc

  implicit none

  type(element_work_t) :: ework
  real(rk), parameter :: young = 200.0e3_rk
  real(rk), parameter :: nu = 0.3_rk

  character(len=512) :: meshes
  character(len=32), parameter :: mesh_names(2) = [character(len=32) :: &
    & 'periodic_box_2x1x1.msh', 'periodic_cube.msh']
  integer :: i_mesh, i_case, unit_log

  ! The mesh data of the current mesh (module-like state of the program).
  type(mesh_t) :: mesh
  type(mesh_file_t) :: ref
  type(element_t) :: elt
  real(rk), allocatable :: k_elt(:, :, :), k_ref_all(:, :, :)
  real(rk) :: lengths(3)

  call par_init()
  call get_command_argument(1, meshes)
  open (newunit=unit_log, status='scratch', action='readwrite', form='formatted')
  call log_set_unit(unit_log)

  do i_mesh = 1, size(mesh_names)
    call load(trim(meshes) // '/' // trim(mesh_names(i_mesh)))
    do i_case = 1, 3
      call check_case(trim(mesh_names(i_mesh)), i_case)
    end do
    call check_stress_average(trim(mesh_names(i_mesh)))
  end do

  call log_set_unit(output_unit)
  close (unit_log)
  call test_summary()
  call par_finalize()

contains

  !> Load a mesh (distributed and whole), the element matrices and the lengths of its drive
  !> pairs (from the whole mesh's coordinates).
  !>
  !> Arguments:
  !>   file  (in)  character  mesh file
  !> Locals:
  !>   mat_s       real(rk)(:, :, :)  [s] at the points [MPa]
  !>   i_elt       integer(ik)        element
  !>   drive_node  integer(ik)(2, 3)  drive pairs
  !>   sigma       real(rk)(3)        their signs
  !>   j           integer            axis
  subroutine load(file)

    character(len=*), intent(in) :: file

    real(rk), allocatable :: mat_s(:, :, :)
    integer(ik) :: i_elt, drive_node(2, 3)
    real(rk) :: sigma(3)
    integer :: j

    call mesh_load_blocks(mesh, file)
    call mesh_file_load(ref, file)
    call element_init(elt, mesh%elt_type)
    allocate (mat_s(6, 6, elt%num_qpts))
    mat_s = spread(isotropic_mandel(young, nu), 3, int(elt%num_qpts))
    if (allocated(k_elt)) deallocate (k_elt, k_ref_all)
    allocate (k_elt(elt%num_dofs, elt%num_dofs, mesh%num_elts))
    do i_elt = 1_ik, mesh%num_elts
      call element_work_alloc(ework, elt)
      call element_stiffness(elt, mesh%coo_s(:, mesh%elt_nodes_loc(:, i_elt)), mat_s, &
        & k_elt(:, :, i_elt), ework)
    end do
    allocate (k_ref_all(elt%num_dofs, elt%num_dofs, ref%num_elts))
    do i_elt = 1_ik, ref%num_elts
      call element_work_alloc(ework, elt)
      call element_stiffness(elt, ref%coo(:, ref%elt_nodes(:, i_elt)), mat_s, &
        & k_ref_all(:, :, i_elt), ework)
    end do
    call periodic_drive_pairs(ref%per_secondary, ref%per_primary, ref%per_vect, drive_node, &
      & sigma)
    do j = 1, 3
      lengths(j) = abs(ref%coo(j, drive_node(1, j)) - ref%coo(j, drive_node(2, j)))
    end do

  end subroutine load

  !> One case of the module header's checks.
  !>
  !> Arguments:
  !>   name    (in)  character  mesh name (messages)
  !>   i_case  (in)  integer    1 small strain (11 22 12), 2 deformation_gradient_rate (11 12 13
  !>                            23), 3 deformation_gradient_rate (11)
  !> Locals:
  !>   label    character(:)        message prefix
  !>   small    logical             small strain
  !>   imposed  logical(9)          imposed components
  !>   rates    real(rk)(9)         imposed rates [1/s]
  !>   cons     type(constraint_t)  constraints (save: gfortran 15 -O3 otherwise warns of
  !>                                uninitialized descriptors, ROADMAP section G)
  !>   num_free integer(ik)         free parameters of Fbar-dot
  !>   num_rows integer(ik)         local rows of the owned dofs, over the ranks
  !>   x, y, tx, tty, ky, d, vel  real(rk)(:)  local vectors
  !>   x_ref, ttkx, e_q, te_q, kte_q, d_ref  real(rk)(:)  whole-mesh vectors
  !>   mask     logical(:)          no prescribed dofs
  !>   coo      real(rk)(:, :)      x = X + v [mm]
  !>   fdot     real(rk)(3, 3)      Fbar-dot from the drive pairs [1/s]
  !>   f_pairs, f_avg  real(rk)(3, 3)  Fbar of x from the drive pairs, and averaged [-]
  !>   jump     real(rk)(3)         a relation's velocity jump minus Fbar-dot (P o l)
  !>   sums, sums_all  real(rk)(3)  sum w v over the owned nodes, over the ranks
  !>   v_ref, w_ref    real(rk)     serial V_0 and int X dOmega_0 (w_ref(3))
  !>   err, err_all, scale, a, b, dq  real(rk)  errors, scales, sides of identities
  !>   i, k, g, q, r, e, s, p  integer(ik)  indices
  !>   ok       logical             all entries agree
  subroutine check_case(name, i_case)

    character(len=*), intent(in) :: name
    integer, intent(in) :: i_case

    character(len=:), allocatable :: label
    logical :: small, imposed(9), ok
    real(rk) :: rates(9), fdot(3, 3), f_pairs(3, 3), f_avg(3, 3), jump(3), sums(3), sums_all(3)
    real(rk) :: v_ref, w_ref(3), err, err_all, scale, a, b, dq
    type(constraint_t), save :: cons
    integer(ik) :: num_free, num_rows, i, k, g, q, r, e, s, p
    real(rk), allocatable :: x(:), y(:), tx(:), tty(:), ky(:), d(:), vel(:)
    real(rk), allocatable :: x_ref(:), ttkx(:), e_q(:), te_q(:), kte_q(:), d_ref(:)
    real(rk), allocatable :: coo(:, :)
    logical, allocatable :: mask(:)

    ! The case: kinematics, imposed components (index 3 (i - 1) + j) and rates.
    label = name // ' case ' // achar(48 + i_case) // ': '
    imposed = .false.
    select case (i_case)
    case (1)
      small = .true.
      imposed([1, 5, 2, 4]) = .true.
    case (2)
      small = .false.
      imposed([1, 2, 3, 6]) = .true.
    case default
      small = .false.
      imposed(1) = .true.
    end select
    rates = 0.0_rk
    do q = 1_ik, 9_ik
      if (imposed(q)) rates(q) = 1.0e-3_rk*sin(1.3_rk*q)
    end do
    if (small) rates(4) = rates(2)
    call periodic_build(cons, mesh, elt, small, imposed, .true.)
    call constraint_set_rates(cons, rates)

    ! Table: one row per secondary dof but one freed dof per free parameter; no constrained
    ! column; the top primary among the columns.
    num_free = int(count(.not. imposed), ik)
    if (small) num_free = int(count(.not. imposed([1, 5, 9, 2, 3, 6])), ik)
    call assert_equal(int(cons%num_rows_glob), int(3_ik*ref%num_periodicity - num_free), &
      & label // 'one row per secondary dof but the freed drive dofs')
    call par_sum(int(count(cons%row_dof <= 3_ik*mesh%num_nodes_owned), ik), num_rows)
    call assert_equal(int(num_rows), int(cons%num_rows_glob), label // 'each row owned once')
    ok = .true.
    do e = 1_ik, int(size(cons%glob_col), ik)
      if (find_sorted(cons%glob_dof, cons%glob_col(e)) > 0_ik) ok = .false.
    end do
    call assert_true(ok, label // 'no column is a constrained dof')
    ok = .true.
    do r = 1_ik, cons%num_rows_glob
      if (.not. any(cons%glob_col(cons%glob_start(r):cons%glob_start(r + 1_ik) - 1_ik) == &
        & cons%glob_prim(r))) ok = .false.
    end do
    call assert_true(ok, label // 'the top primary is a column')

    ! Relations: v = T v_hat + g from a reduced vector of the global ids.
    allocate (vel(3*mesh%num_nodes_loc), mask(3*mesh%num_nodes_loc))
    mask = .false.
    do i = 1_ik, mesh%num_nodes_loc
      do k = 1_ik, 3_ik
        vel(3*(i - 1) + k) = 1.0e-3_rk*sin(0.1_rk*mesh%node_id_glob(i) + k)
      end do
    end do
    where (cons%is_constrained) vel = 0.0_rk
    call constraint_expand(cons, vel, .true.)
    fdot = constraint_macro_rate(cons, vel)
    ! (Fbar-dot row by row is the index order 3 (i - 1) + j.)
    err = maxval(abs(reshape(transpose(fdot), [9]) - rates), mask=imposed)
    call assert_true(err <= 1e-12_rk*maxval(abs(rates)), label // 'Fbar-dot = the imposed ' // &
      & 'rates on the imposed components')
    if (small) call assert_true(maxval(abs(fdot - transpose(fdot))) <= 1e-12_rk* &
      & maxval(abs(fdot)), label // 'small strain: Fbar-dot symmetric')
    err = 0.0_rk
    scale = 0.0_rk
    do q = 1_ik, ref%num_periodicity
      s = node_loc_from_glob(mesh, ref%per_secondary(q))
      p = node_loc_from_glob(mesh, ref%per_primary(q))
      if (s == 0_ik .or. p == 0_ik) cycle
      jump = vel(3*s - 2:3*s) - vel(3*p - 2:3*p)
      scale = max(scale, maxval(abs(jump)))
      jump = jump - matmul(fdot, real(ref%per_vect(:, q), rk)*lengths)
      err = max(err, maxval(abs(jump)))
    end do
    call par_max(err, err_all)
    call par_max(scale, a)
    call assert_true(err_all <= 1e-12_rk*a, label // 'relations v(s) - v(p) = Fbar-dot (P o l)')

    ! V8: Fbar of x = X + v from the drive pairs and as the volume average of F.
    coo = mesh%coo_s + reshape(vel, [3_ik, mesh%num_nodes_loc])
    f_pairs = constraint_macro_rate(cons, reshape(coo, [3_ik*mesh%num_nodes_loc]))
    f_avg = periodic_fbar_average(mesh, elt, coo)
    call assert_true(maxval(abs(f_pairs - f_avg)) <= 1e-10_rk*maxval(abs(fdot)), label // &
      & 'Fbar from the drive pairs = volume average of F')

    ! Weights: V_0 and W against a serial quadrature.
    call serial_moments(v_ref, w_ref)
    call assert_close(cons%volume0, v_ref, label // 'sum of the weights = V_0', rel_tol=1e-12_rk)
    call assert_true(maxval(abs(cons%moment0 - w_ref)) <= 1e-12_rk*maxval(abs(w_ref)), &
      & label // 'sum w X = int X dOmega_0')

    ! Zero-mean shift: sum w v = Fbar-dot W afterwards.
    call constraint_zero_mean(cons, mesh, vel)
    sums = 0.0_rk
    do i = 1_ik, mesh%num_nodes_owned
      sums = sums + cons%node_weight(i)*vel(3*i - 2:3*i)
    end do
    call par_sum(sums, sums_all)
    call assert_true(maxval(abs(sums_all - matmul(fdot, cons%moment0))) <= 1e-12_rk* &
      & maxval(abs(fdot))*maxval(abs(cons%moment0)) + 1e-12_rk*maxval(abs(sums_all)), &
      & label // 'zero-mean fluctuation: sum w v = Fbar-dot W')

    ! Adjointness: x reduced, y arbitrary.
    allocate (x(3*mesh%num_nodes_loc), y(3*mesh%num_nodes_loc))
    do i = 1_ik, mesh%num_nodes_loc
      do k = 1_ik, 3_ik
        x(3*(i - 1) + k) = sin(0.1_rk*mesh%node_id_glob(i) + k)
        y(3*(i - 1) + k) = cos(0.07_rk*mesh%node_id_glob(i) + 2*k)
      end do
    end do
    where (cons%is_constrained) x = 0.0_rk
    tx = x
    call constraint_expand(cons, tx, .false.)
    tty = y
    call constraint_reduce(cons, 3_ik*mesh%num_nodes_owned, tty)
    a = dist_dot(tx, y, 3_ik*mesh%num_nodes_owned)
    b = dist_dot(x, tty, 3_ik*mesh%num_nodes_owned)
    call assert_close(a, b, label // 'T and T^T adjoint', rel_tol=1e-12_rk)

    ! Reduced operator against the serial T^T K T x.
    allocate (x_ref(3*ref%num_nodes))
    do g = 1_ik, ref%num_nodes
      x_ref(3*g - 2:3*g) = [(sin(0.1_rk*g + k), k = 1_ik, 3_ik)]
    end do
    x_ref(cons%glob_dof) = 0.0_rk
    ttkx = serial_ttkt(cons, x_ref)
    call halo_forward(mesh, x, 3_ik)
    allocate (ky(size(x)))
    call constraint_matvec(cons, mesh, k_elt, mask, x, ky)
    scale = maxval(abs(ttkx))
    ok = .true.
    do i = 1_ik, mesh%num_nodes_owned
      g = mesh%node_id_glob(i)
      if (any(abs(ky(3*i - 2:3*i) - ttkx(3*g - 2:3*g)) > 1e-12_rk*scale)) ok = .false.
    end do
    call assert_true(ok, label // 'reduced operator = serial T^T K T')

    ! Exact diagonal: every column of the table and every 23rd dof.
    allocate (d(3*mesh%num_nodes_loc), e_q(3*ref%num_nodes))
    call constraint_diagonal(cons, mesh, k_elt, d)
    err = 0.0_rk
    scale = 0.0_rk
    do q = 1_ik, 3_ik*ref%num_nodes
      if (find_sorted(cons%glob_dof, q) > 0_ik) cycle
      if (find_sorted(cons%all_cols, q) == 0_ik .and. mod(q, 23_ik) /= 0_ik) cycle
      e_q = 0.0_rk
      e_q(q) = 1.0_rk
      te_q = serial_expand(cons, e_q)
      kte_q = serial_k(te_q)
      dq = dot_product(te_q, kte_q)
      scale = max(scale, abs(dq))
      i = owned_local((q - 1_ik)/3_ik + 1_ik)
      if (i > 0_ik) err = max(err, abs(d(3*(i - 1_ik) + mod(q - 1_ik, 3_ik) + 1_ik) - dq))
    end do
    call par_max(err, a)
    call assert_true(a <= 1e-12_rk*scale, label // 'exact diagonal = (T e_q)^T K (T e_q)')

    ! The lumped diagonal: the element diagonals onto the top primaries.
    call constraint_diagonal_lumped(cons, mesh, k_elt, d)
    d_ref = serial_diagonal_lumped(cons)
    err = 0.0_rk
    do i = 1_ik, mesh%num_nodes_owned
      g = mesh%node_id_glob(i)
      err = max(err, maxval(abs(d(3*i - 2:3*i) - d_ref(3*g - 2:3*g))))
    end do
    call par_max(err, a)
    call assert_true(a <= 1e-12_rk*maxval(abs(d_ref)), label // 'lumped diagonal = serial ' // &
      & '(T+)^T diag K')

  end subroutine check_case

  !> V8: the volume average of a synthetic stress field (periodic_stress_average) against a serial
  !> quadrature on the whole mesh.
  !>
  !> Arguments:
  !>   name  (in)  character  mesh name (messages)
  !> Locals:
  !>   stress_s  real(rk)(:, :, :)  (6, num_qpts, num_elts) Mandel stress of the local elements
  !>   avg, ref_avg  real(rk)(6)   the average (tensor components), the serial one
  !>   sums      real(rk)(7)       serial sums of w det J sigma and of the volume
  !>   dshape_s  real(rk)(:, :)    shape derivatives (unused)
  !>   det_jac   real(rk)          det J
  !>   i_elt, i_qpt  integer(ik)   element, point
  subroutine check_stress_average(name)

    character(len=*), intent(in) :: name

    real(rk), allocatable :: stress_s(:, :, :), dshape_s(:, :)
    real(rk) :: avg(6), ref_avg(6), sums(7), det_jac
    integer(ik) :: i_elt, i_qpt

    allocate (stress_s(6, elt%num_qpts, mesh%num_elts), dshape_s(3, elt%num_nodes))
    do i_elt = 1_ik, mesh%num_elts
      do i_qpt = 1_ik, elt%num_qpts
        stress_s(:, i_qpt, i_elt) = synthetic_stress(mesh%elt_id_glob(i_elt), i_qpt)
      end do
    end do
    avg = periodic_stress_average(mesh, elt, mesh%coo_s, stress_s)
    sums = 0.0_rk
    do i_elt = 1_ik, ref%num_elts
      do i_qpt = 1_ik, elt%num_qpts
        call element_geometry_qpt(elt, ref%coo(:, ref%elt_nodes(:, i_elt)), i_qpt, dshape_s, &
          & det_jac)
        sums = sums + elt%qpt_weight(i_qpt)*det_jac*[synthetic_stress(i_elt, i_qpt), 1.0_rk]
      end do
    end do
    ref_avg = sums(1:6)/sums(7)
    ref_avg(4:6) = ref_avg(4:6)/sqrt2
    call assert_true(maxval(abs(avg - ref_avg)) <= 1e-12_rk*maxval(abs(ref_avg)), name // &
      & ': volume-averaged stress = serial quadrature')

  end subroutine check_stress_average

  !> A synthetic Mandel stress of an element and point [MPa].
  !>
  !> Arguments:
  !>   id     (in)  integer(ik)  global element id
  !>   i_qpt  (in)  integer(ik)  point
  !>   s      (out, result) real(rk)(6)
  !> Locals:
  !>   c  integer  component
  pure function synthetic_stress(id, i_qpt) result(s)

    integer(ik), intent(in) :: id, i_qpt
    real(rk) :: s(6)

    integer :: c

    s = [(100.0_rk*sin(0.3_rk*id + 0.7_rk*i_qpt + c), c = 1, 6)]

  end function synthetic_stress

  !> Serial V_0 and int X dOmega_0 of the whole mesh by the elements' quadrature.
  !>
  !> Arguments:
  !>   v_ref  (out)  real(rk)     volume [mm^3]
  !>   w_ref  (out)  real(rk)(3)  int X dOmega_0 [mm^4]
  !> Locals:
  !>   dshape_s  real(rk)(:, :)  shape derivatives (unused)
  !>   det_jac   real(rk)        det J
  !>   i_elt, i_qpt  integer(ik) element, point
  subroutine serial_moments(v_ref, w_ref)

    real(rk), intent(out) :: v_ref, w_ref(3)

    real(rk), allocatable :: dshape_s(:, :)
    real(rk) :: det_jac
    integer(ik) :: i_elt, i_qpt

    allocate (dshape_s(3, elt%num_nodes))
    v_ref = 0.0_rk
    w_ref = 0.0_rk
    do i_elt = 1_ik, ref%num_elts
      do i_qpt = 1_ik, elt%num_qpts
        call element_geometry_qpt(elt, ref%coo(:, ref%elt_nodes(:, i_elt)), i_qpt, dshape_s, &
          & det_jac)
        v_ref = v_ref + elt%qpt_weight(i_qpt)*det_jac
        w_ref = w_ref + elt%qpt_weight(i_qpt)*det_jac*matmul(ref%coo(:, ref%elt_nodes(:, &
          & i_elt)), elt%shape_qpt(:, i_qpt))
      end do
    end do

  end subroutine serial_moments

  !> Serial lumped periodic diagonal: every element diagonal K_ii onto its dof, or onto the top
  !> primary of its row when constrained; zero at the constrained dofs.
  !>
  !> Arguments:
  !>   cons  (in)  type(constraint_t)  table
  !>   d     (out, result) real(rk)(:)  (3 num_nodes)
  !> Locals:
  !>   e, a, c  integer      element, node, component
  !>   dof, r   integer(ik)  global dof, its row
  function serial_diagonal_lumped(cons) result(d)

    type(constraint_t), intent(in) :: cons
    real(rk), allocatable :: d(:)

    integer :: e, a, c
    integer(ik) :: dof, r

    allocate (d(3*ref%num_nodes))
    d = 0.0_rk
    do e = 1, int(ref%num_elts)
      do a = 1, int(elt%num_nodes)
        do c = 1, 3
          dof = 3_ik*(ref%elt_nodes(a, e) - 1_ik) + c
          r = find_sorted(cons%glob_dof, dof)
          if (r > 0_ik) dof = cons%glob_prim(r)
          d(dof) = d(dof) + k_ref_all(3*(a - 1) + c, 3*(a - 1) + c, e)
        end do
      end do
    end do
    d(cons%glob_dof) = 0.0_rk

  end function serial_diagonal_lumped

  !> Serial T x (no offsets) on the whole mesh, from the replicated table.
  !>
  !> Arguments:
  !>   cons  (in)  type(constraint_t)  table
  !>   x     (in)  real(rk)(:)         (3 num_nodes) independent values (zero elsewhere)
  !>   tx    (out, result) real(rk)(:)
  !> Locals:
  !>   r, e  integer(ik)  row, column entry
  function serial_expand(cons, x) result(tx)

    type(constraint_t), intent(in) :: cons
    real(rk), intent(in) :: x(:)
    real(rk), allocatable :: tx(:)

    integer(ik) :: r, e

    tx = x
    do r = 1_ik, cons%num_rows_glob
      tx(cons%glob_dof(r)) = 0.0_rk
      do e = cons%glob_start(r), cons%glob_start(r + 1_ik) - 1_ik
        tx(cons%glob_dof(r)) = tx(cons%glob_dof(r)) + cons%glob_coef(e)*x(cons%glob_col(e))
      end do
    end do

  end function serial_expand

  !> Serial K x on the whole mesh (element by element, k_ref_all).
  !>
  !> Arguments:
  !>   x   (in)  real(rk)(:)  (3 num_nodes) vector
  !>   kx  (out, result) real(rk)(:)
  !> Locals:
  !>   e, a   integer          element, node
  !>   dofs   integer(ik)(:)   element dofs
  function serial_k(x) result(kx)

    real(rk), intent(in) :: x(:)
    real(rk), allocatable :: kx(:)

    integer :: e, a
    integer(ik), allocatable :: dofs(:)

    allocate (kx(size(x)), dofs(elt%num_dofs))
    kx = 0.0_rk
    do e = 1, int(ref%num_elts)
      do a = 1, int(elt%num_nodes)
        dofs(3*a - 2:3*a) = 3_ik*(ref%elt_nodes(a, e) - 1_ik) + [1_ik, 2_ik, 3_ik]
      end do
      kx(dofs) = kx(dofs) + matmul(k_ref_all(:, :, e), x(dofs))
    end do

  end function serial_k

  !> Serial T^T K T x on the whole mesh.
  !>
  !> Arguments:
  !>   cons  (in)  type(constraint_t)  table
  !>   x     (in)  real(rk)(:)         (3 num_nodes) independent values
  !>   y     (out, result) real(rk)(:)
  !> Locals:
  !>   r, e  integer(ik)  row, column entry
  function serial_ttkt(cons, x) result(y)

    type(constraint_t), intent(in) :: cons
    real(rk), intent(in) :: x(:)
    real(rk), allocatable :: y(:)

    integer(ik) :: r, e

    y = serial_k(serial_expand(cons, x))
    do r = 1_ik, cons%num_rows_glob
      do e = cons%glob_start(r), cons%glob_start(r + 1_ik) - 1_ik
        y(cons%glob_col(e)) = y(cons%glob_col(e)) + cons%glob_coef(e)*y(cons%glob_dof(r))
      end do
      y(cons%glob_dof(r)) = 0.0_rk
    end do

  end function serial_ttkt

  !> Local number of a node owned by this rank (0 otherwise).
  !>
  !> Arguments:
  !>   id     (in)  integer(ik)  global node id
  !>   i_loc  (out, result) integer(ik)
  function owned_local(id) result(i_loc)

    integer(ik), intent(in) :: id
    integer(ik) :: i_loc

    i_loc = node_loc_from_glob(mesh, id)
    if (i_loc > mesh%num_nodes_owned) i_loc = 0_ik

  end function owned_local

  !> Isotropic stiffness, Mandel [MPa].
  !>
  !> Arguments:
  !>   e_mod, poisson  (in)  real(rk)  Young's modulus [MPa], Poisson's ratio
  !>   stiff  (out, result) real(rk)(6, 6)
  !> Locals:
  !>   lam, mu  real(rk)  Lame constants [MPa]
  !>   i, j     integer   indices
  pure function isotropic_mandel(e_mod, poisson) result(stiff)

    real(rk), intent(in) :: e_mod, poisson
    real(rk) :: stiff(6, 6)

    real(rk) :: lam, mu
    integer :: i, j

    lam = e_mod*poisson/((1.0_rk + poisson)*(1.0_rk - 2.0_rk*poisson))
    mu = e_mod/(2.0_rk*(1.0_rk + poisson))
    stiff = 0.0_rk
    do i = 1, 3
      do j = 1, 3
        stiff(i, j) = lam
      end do
    end do
    do i = 1, 6
      stiff(i, i) = stiff(i, i) + 2.0_rk*mu
    end do

  end function isotropic_mandel

end program test_periodic
