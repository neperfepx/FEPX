! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> MPI test of sim_writer_mod (1, 2 and 4 ranks): writes every variable of the step 4.3 subset
!> for two steps on the distributed sx2 mesh, then reads every file back on the root and checks
!> it against values computed from global ids and the whole mesh, so the files must be identical
!> whatever the number of ranks:
!>   - node files: one line per node in global id order, 3 values, exact (es24.16 round trip);
!>   - element files: one line per element in global id order, the 6 tensor components (Mandel
!>     shear divided by sqrt 2), exact;
!>   - force files: v2's two header lines (identical text), then one row per call;
!>   - no .sim summary (D109); the copied inputs, byte for byte;
!>   - variables not printed are not written;
!>   - values are placed by global id: disp and strain_el are written through a copy of the
!>     mesh whose owned nodes and elements are listed in reverse order;
!>   - slip-system variables (step 6.5, sim_write_elts_slip): crss with 2 or 3 values per
!>     element (by the parity of the global id, as for two phases), written through the reversed
!>     copy: each line has its element's count of values, exact;
!>   - variables of k values as given (sim_write_elts_values): work, one value per line, exact;
!>   - the convergence record (step 6.5b, sim_write_convergence): v2's header, then one row per
!>     global iteration of each call (increment, iteration, method, six norms, CG iterations),
!>     exact; listed among the element results of the .sim file, with no results/elts
!>     directory;
!>   - the sub-incrementation columns (Q35, D63; sim_open attempts): the header with
!>     "attempt A|R", and every row with its attempt and A or R, exact;
!>   - the boundary-velocity iteration columns (step 8.2b; sim_open bc_iterations): the header
!>     with "attempt A|R  bc_iter    solve", and every row with its attempt, A, iteration and
!>     solve (0 or the perturbed component), exact.
!>
!> Source : FEPX (ROADMAP step 4.3). The first argument is the mesh directory
!>          (tests/regression/cases/meshes).
!>
!> Locals:
!>   meshes       character(512)        mesh directory
!>   mesh         type(mesh_t)          rank-local mesh
!>   mesh_rev     type(mesh_t)          the same, owned nodes and elements in reverse order
!>   n_own        integer(ik)           owned nodes
!>   ref          type(mesh_file_t)     whole mesh (root), for the expected coordinates
!>   sim          type(sim_writer_t)    the writer
!>   dir, cfg     character(:), allocatable  output directory and input file of this run
!>   inputs       character(512)(2)     input files given to the writer
!>   unit_log     integer               scratch log unit
!>   unit, ios    integer               file unit, status
!>   step         integer               step index
!>   i, j         integer(ik)           indices
!>   disp         real(rk)(:, :)        (3, num_nodes_loc) node values
!>   stress       real(rk)(:, :)        (6, num_elts) Mandel element values
!>   slipv        real(rk)(:, :)        (3, num_elts) slip-system element values
!>   counts       integer(ik)(:)        (num_elts) slip systems per element
!>   workv        real(rk)(:, :)        (1, num_elts) scalar element values
!>   loads        real(rk)(:, :)        (3, num_fasets) face-set forces
!>   areas        real(rk)(:)           (num_fasets) face-set areas
!>   line         character(1024)       a line read back
!>   row          real(rk)(6)           values read back
!>   expect       real(rk)(6)           expected element value
!>   step_text    character(8)          step number as text
!>   ok           logical               all values of a file as expected
!>   ints         integer(2)            step and increment read back
!>   exists       logical               whether a file exists
!>   norms        real(rk)(6, 2)        norms of two global iterations
!>   method       character(2)          method read back
!>   methods      character(2)(2)       methods of two global iterations
!>   cgs          integer(2)            CG iterations read back (and the iteration number)
!>   attempt_read integer               attempt read back (sub-incrementation columns)
!>   bc_read      integer(2)            iteration and solve read back (boundary-velocity
!>                                      iteration columns)
!>   mark         character(1)          A or R read back
!>   expect_rows  integer(3, 3)         expected iteration, CG count index and attempt per row
program test_sim_writer

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk, ik
  use constants_mod, only: sqrt2
  use log_mod, only: log_set_unit
  use parallel_mod, only: par_init, par_finalize, par_barrier, num_ranks, my_rank, root_rank
  use mesh_file_mod, only: mesh_file_t, mesh_file_load
  use mesh_mod, only: mesh_t
  use mesh_blocks_mod, only: mesh_load_blocks
  use sim_writer_mod, only: sim_writer_t, sim_open, sim_write_nodes, sim_write_elts, &
    & sim_write_elts_values, sim_write_elts_slip, sim_write_forces, sim_write_convergence, &
    & sim_close
  use assert_mod, only: assert_true, assert_equal, test_summary

  implicit none

  character(len=512) :: meshes
  type(mesh_t) :: mesh, mesh_rev
  integer(ik) :: n_own
  type(mesh_file_t) :: ref
  type(sim_writer_t) :: sim
  character(len=:), allocatable :: dir, cfg
  character(len=512) :: inputs(2)
  integer :: unit_log, unit, ios, step
  integer(ik) :: i, j
  real(rk), allocatable :: disp(:, :), stress(:, :), loads(:, :), areas(:), slipv(:, :), &
    & workv(:, :)
  integer(ik), allocatable :: counts(:)
  character(len=1024) :: line
  real(rk) :: row(6), expect(6)
  character(len=8) :: step_text
  logical :: ok, exists
  integer :: ints(2), cgs(2)
  real(rk) :: norms(6, 2)
  character(len=2) :: method, methods(2)
  integer :: attempt_read, expect_rows(3, 3), bc_read(2)
  character(len=1) :: mark

  call par_init()
  call get_command_argument(1, meshes)

  ! Mesh (log to a scratch unit); the whole mesh on the root for the expected coordinates.
  open (newunit=unit_log, status='scratch', action='readwrite', form='formatted')
  call log_set_unit(unit_log)
  call mesh_load_blocks(mesh, trim(meshes) // '/sx2.msh')
  if (my_rank == root_rank) call mesh_file_load(ref, trim(meshes) // '/sx2.msh')
  call log_set_unit(output_unit)
  close (unit_log)

  ! Fresh output directory and a small input file for this rank count.
  write (step_text, '(i0)') num_ranks
  dir = 'test_sim_np' // trim(step_text) // '.sim'
  cfg = 'test_sim_np' // trim(step_text) // '.cfg'
  if (my_rank == root_rank) then
    call execute_command_line('rm -rf ' // dir)
    open (newunit=unit, file=cfg, status='replace', action='write')
    write (unit, '(a)') '# input copied by the writer'
    write (unit, '(a)') 'print coo'
    close (unit)
  end if
  call par_barrier()

  ! Write: every variable but disp at step 0 (disp written for step 1 only), two force rows.
  ! (The input list is a fixed-length variable: gfortran 15 builds a typed array constructor
  ! whose first item is a deferred-length string with that item's length; ROADMAP section G.)
  inputs(1) = cfg
  inputs(2) = trim(meshes) // '/sx2.msh'
  call sim_open(sim, dir, [character(len=32) :: 'coo', 'stress', 'disp', 'forces', &
    & 'strain_el', 'crss', 'work', 'convergence'], mesh, inputs)
  allocate (disp(3, mesh%num_nodes_loc), stress(6, mesh%num_elts))
  allocate (slipv(3, mesh%num_elts), counts(mesh%num_elts), workv(1, mesh%num_elts))

  ! The reversed copy (the writer must place values by global id, not by gathered position).
  mesh_rev = mesh
  n_own = mesh%num_nodes_owned
  mesh_rev%node_id_glob(1:n_own) = mesh%node_id_glob(n_own:1:-1)
  mesh_rev%elt_id_glob = mesh%elt_id_glob(mesh%num_elts:1:-1)
  allocate (loads(3, size(mesh%fasets)), areas(size(mesh%fasets)))
  do step = 0, 1
    ! Node values from the reversed copy's global ids; ghost columns get a value that must not
    ! be written.
    do i = 1_ik, mesh%num_nodes_loc
      disp(:, i) = node_value(mesh_rev%node_id_glob(i), step)
      if (i > n_own) disp(:, i) = -1.0e300_rk
    end do
    call sim_write_nodes(sim, mesh, 'coo', mesh%coo_s, step)
    if (step == 1) call sim_write_nodes(sim, mesh_rev, 'disp', disp, step)
    call sim_write_nodes(sim, mesh, 'vel', disp, step)

    ! Element values from global ids: stress through the mesh, strain_el through the copy.
    do i = 1_ik, mesh%num_elts
      stress(:, i) = elt_value(mesh%elt_id_glob(i), step)
    end do
    call sim_write_elts(sim, mesh, 'stress', stress, step)
    call sim_write_elts(sim, mesh, 'strain', stress, step)
    do i = 1_ik, mesh%num_elts
      stress(:, i) = elt_value(mesh_rev%elt_id_glob(i), step)
    end do
    call sim_write_elts(sim, mesh_rev, 'strain_el', -stress, step)

    ! Slip-system values through the copy (2 or 3 per element), scalars through the mesh.
    do i = 1_ik, mesh%num_elts
      expect = elt_value(mesh_rev%elt_id_glob(i), step)
      slipv(:, i) = expect(1:3)
      counts(i) = 2_ik + mod(mesh_rev%elt_id_glob(i), 2_ik)
      expect = elt_value(mesh%elt_id_glob(i), step + 1)
      workv(1, i) = expect(1)
    end do
    call sim_write_elts_slip(sim, mesh_rev, 'crss', slipv, counts, step)
    call sim_write_elts_values(sim, mesh, 'work', workv, step)
    do i = 1_ik, size(mesh%fasets, kind=ik)
      loads(:, i) = [1.0_rk, -2.0_rk, 3.0_rk]*(i + step/7.0_rk)
      areas(i) = 1.0_rk + i/3.0_rk
    end do
    call sim_write_forces(sim, step, 5*step, loads, areas, 0.25_rk*step)

    ! Convergence record: step + 1 global iterations of increment 3 step (SA, then NR).
    do i = 1_ik, 2_ik
      norms(:, i) = [(real(j, rk)/(7.0_rk*i + step), j = 1_ik, 6_ik)]
    end do
    methods = [character(len=2) :: 'SA', 'NR']
    cgs = [11, 22]
    call sim_write_convergence(sim, 3*step, methods(1:step + 1), norms(:, 1:step + 1), &
      & cgs(1:step + 1))
  end do
  call sim_close(sim)
  call par_barrier()

  ! Writer state on every rank: the variables by entity (forces with the elements, as v2), the
  ! completed steps, and force files open on the root only.
  call assert_equal(size(sim%node_vars), 2, 'node variables')
  call assert_equal(size(sim%elt_vars), 6, 'element variables (forces, convergence included)')
  call assert_equal(sim%num_steps, 1, 'completed steps')
  call assert_equal(size(sim%force_units), merge(size(mesh%fasets), 0, my_rank == root_rank), &
    & 'force files on the root only')

  ! --- Read back on the root ---
  if (my_rank == root_rank) then
    ! Coordinates (both steps), displacement (step 1 only), not-printed variables absent.
    do step = 0, 1
      write (step_text, '(i0)') step
      open (newunit=unit, file=dir // '/results/nodes/coo/coo.step' // trim(step_text), &
        & status='old', action='read')
      ok = .true.
      do i = 1_ik, ref%num_nodes
        read (unit, *) row(1:3)
        ok = ok .and. all(abs(row(1:3) - ref%coo(:, i)) <= 0.0_rk)
      end do
      read (unit, '(a)', iostat=ios) line
      call assert_true(ok .and. ios /= 0, 'coo step ' // trim(step_text) // ': all nodes, exact')
      close (unit)
    end do
    inquire (file=dir // '/results/nodes/disp/disp.step0', exist=exists)
    call assert_true(.not. exists, 'disp step 0 not written')
    open (newunit=unit, file=dir // '/results/nodes/disp/disp.step1', status='old', &
      & action='read')
    ok = .true.
    do i = 1_ik, ref%num_nodes
      read (unit, '(a)') line
      read (line, *) row(1:3)
      ok = ok .and. all(abs(row(1:3) - node_value(i, 1)) <= 0.0_rk)
      ok = ok .and. count_values(line) == 3
    end do
    call assert_true(ok, 'disp step 1: global id order, 3 values per line, exact')
    close (unit)
    inquire (file=dir // '/results/nodes/vel', exist=exists)
    call assert_true(.not. exists, 'vel (not printed) not written')
    inquire (file=dir // '/results/elts/strain', exist=exists)
    call assert_true(.not. exists, 'strain (not printed) not written')

    ! Stress and elastic strain: tensor components, shear = Mandel / sqrt 2.
    open (newunit=unit, file=dir // '/results/elts/stress/stress.step1', status='old', &
      & action='read')
    ok = .true.
    do i = 1_ik, ref%num_elts
      read (unit, '(a)') line
      read (line, *) row
      ok = ok .and. count_values(line) == 6
      expect = elt_value(i, 1)
      ok = ok .and. all(abs(row(1:3) - expect(1:3)) <= 0.0_rk)
      ok = ok .and. all(abs(row(4:6) - elt_value_shear(i, 1)) <= 0.0_rk)
    end do
    read (unit, '(a)', iostat=ios) line
    call assert_true(ok .and. ios /= 0, 'stress step 1: all elements, tensor components, exact')
    close (unit)
    open (newunit=unit, file=dir // '/results/elts/strain_el/strain_el.step0', status='old', &
      & action='read')
    ok = .true.
    do i = 1_ik, ref%num_elts
      read (unit, *) row
      expect = elt_value(i, 0)
      ok = ok .and. all(abs(row(1:3) + expect(1:3)) <= 0.0_rk)
      ok = ok .and. all(abs(row(4:6) + elt_value_shear(i, 0)) <= 0.0_rk)
    end do
    call assert_true(ok, 'strain_el step 0: exact')
    close (unit)

    ! Slip-system values: each element's count, by global id, exact; scalars, one per line.
    open (newunit=unit, file=dir // '/results/elts/crss/crss.step1', status='old', &
      & action='read')
    ok = .true.
    do i = 1_ik, ref%num_elts
      read (unit, '(a)') line
      j = 2_ik + mod(i, 2_ik)
      ok = ok .and. count_values(line) == j
      read (line, *) row(1:j)
      expect = elt_value(i, 1)
      ok = ok .and. all(abs(row(1:j) - expect(1:j)) <= 0.0_rk)
    end do
    read (unit, '(a)', iostat=ios) line
    call assert_true(ok .and. ios /= 0, 'crss step 1: each element''s count, global id, exact')
    close (unit)
    open (newunit=unit, file=dir // '/results/elts/work/work.step0', status='old', &
      & action='read')
    ok = .true.
    do i = 1_ik, ref%num_elts
      read (unit, '(a)') line
      ok = ok .and. count_values(line) == 1
      read (line, *) row(1)
      expect = elt_value(i, 1)
      ok = ok .and. abs(row(1) - expect(1)) <= 0.0_rk
    end do
    call assert_true(ok, 'work step 0: one value per element, exact')
    close (unit)

    ! Force files: v2's header, then step, incr, force, area, time.
    do j = 1_ik, size(ref%fasets, kind=ik)
      open (newunit=unit, file=dir // '/results/forces/' // ref%fasets(j)%label, &
        & status='old', action='read')
      read (unit, '(a)') line
      call assert_equal(trim(line), '% step     incr     Fx             Fy             ' // &
        & 'Fz             area a         time', 'force header 1 ' // ref%fasets(j)%label)
      read (unit, '(a)') line
      call assert_equal(trim(line), '% -------- -------- -------------- -------------- ' // &
        & '-------------- -------------- --------------', 'force header 2')
      ok = .true.
      do step = 0, 1
        read (unit, *) ints, row(1:5)
        ok = ok .and. ints(1) == step .and. ints(2) == 5*step
        ok = ok .and. all(abs(row(1:3) - [1.0_rk, -2.0_rk, 3.0_rk]*(j + step/7.0_rk)) <= 0.0_rk)
        ok = ok .and. abs(row(4) - (1.0_rk + j/3.0_rk)) <= 0.0_rk
        ok = ok .and. abs(row(5) - 0.25_rk*step) <= 0.0_rk
      end do
      read (unit, '(a)', iostat=ios) line
      call assert_true(ok .and. ios /= 0, 'force rows ' // ref%fasets(j)%label)
      close (unit)
    end do

    ! Convergence file: v2's header, then increment, iteration, method, norms, CG iterations;
    ! no element-result directory.
    open (newunit=unit, file=dir // '/results/convergence', status='old', action='read')
    read (unit, '(a)') line
    call assert_equal(trim(line), '%   incr     iter       NR     r_norm        rx_norm' // &
      & '       f_norm        delu_norm     delux_norm    u_norm      cg_iter', &
      & 'convergence header')
    ok = .true.
    do step = 0, 1
      do i = 1_ik, int(step + 1, ik)
        read (unit, *) ints, method, row(1:6), cgs(1)
        ok = ok .and. ints(1) == 3*step .and. ints(2) == i
        ok = ok .and. method == merge('SA', 'NR', i == 1_ik) .and. cgs(1) == 11*i
        ok = ok .and. all(abs(row(1:6) - [(real(j, rk)/(7.0_rk*i + step), j = 1_ik, 6_ik)]) &
          & <= 0.0_rk)
      end do
    end do
    read (unit, '(a)', iostat=ios) line
    call assert_true(ok .and. ios /= 0, 'convergence rows: one per iteration, exact')
    close (unit)
    inquire (file=dir // '/results/elts/convergence', exist=exists)
    call assert_true(.not. exists, 'no results/elts/convergence directory')

    ! No summary file (D109).
    inquire (file=dir // '/.sim', exist=exists)
    call assert_true(.not. exists, 'no .sim summary')

    ! Copied inputs, byte for byte.
    call assert_true(same_bytes(cfg, dir // '/inputs/' // cfg), 'cfg copied')
    call assert_true(same_bytes(trim(meshes) // '/sx2.msh', dir // '/inputs/sx2.msh'), &
      & 'msh copied')
  end if

  ! Sub-incrementation columns (Q35, D63): a writer with attempts, increment 4 with a rejected
  ! attempt 1 (SA, NR) and an accepted attempt 2 (SA); rows read back exactly.
  write (step_text, '(i0)') num_ranks
  dir = 'test_sim_attempts_np' // trim(step_text) // '.sim'
  if (my_rank == root_rank) call execute_command_line('rm -rf ' // dir)
  call par_barrier()
  call sim_open(sim, dir, [character(len=32) :: 'convergence'], mesh, inputs, attempts=.true.)
  call sim_write_convergence(sim, 4, methods, norms, cgs, attempt=1, accepted=.false.)
  call sim_write_convergence(sim, 4, methods(1:1), norms(:, 1:1), cgs(1:1), attempt=2, &
    & accepted=.true.)
  call sim_close(sim)
  call par_barrier()
  if (my_rank == root_rank) then
    open (newunit=unit, file=dir // '/results/convergence', status='old', action='read')
    read (unit, '(a)') line
    call assert_equal(trim(line), '%   incr     iter       NR     r_norm        rx_norm' // &
      & '       f_norm        delu_norm     delux_norm    u_norm      cg_iter  attempt A|R', &
      & 'convergence header with the attempt columns')
    ! Rows: (iteration, method/norm/CG index, attempt).
    expect_rows = reshape([1, 1, 1, 2, 2, 1, 1, 1, 2], [3, 3])
    ok = .true.
    do i = 1_ik, 3_ik
      read (unit, *) ints, method, row(1:6), cgs(1), attempt_read, mark
      j = int(expect_rows(2, i), ik)
      ok = ok .and. ints(1) == 4 .and. ints(2) == expect_rows(1, i)
      ok = ok .and. method == methods(j) .and. all(abs(row(1:6) - norms(:, j)) <= 0.0_rk)
      ok = ok .and. attempt_read == expect_rows(3, i)
      ok = ok .and. mark == merge('R', 'A', expect_rows(3, i) == 1)
    end do
    read (unit, '(a)', iostat=ios) line
    call assert_true(ok .and. ios /= 0, 'convergence rows with attempt and A|R, exact')
    close (unit)
  end if

  ! Boundary-velocity iteration columns (step 8.2b): increment 5, attempt 1, iteration 1 at the
  ! face velocities (solve 0, SA and NR rows) and with face 22 perturbed (solve 22, SA row).
  dir = 'test_sim_bc_np' // trim(step_text) // '.sim'
  if (my_rank == root_rank) call execute_command_line('rm -rf ' // dir)
  call par_barrier()
  call sim_open(sim, dir, [character(len=32) :: 'convergence'], mesh, inputs, &
    & bc_iterations=.true.)
  call sim_write_convergence(sim, 5, methods, norms, cgs, attempt=1, accepted=.true., bc_iter=1, &
    & bc_solve=0)
  call sim_write_convergence(sim, 5, methods(1:1), norms(:, 1:1), cgs(1:1), attempt=1, &
    & accepted=.true., bc_iter=1, bc_solve=22)
  call sim_close(sim)
  call par_barrier()
  if (my_rank == root_rank) then
    open (newunit=unit, file=dir // '/results/convergence', status='old', action='read')
    read (unit, '(a)') line
    call assert_equal(trim(line), '%   incr     iter       NR     r_norm        rx_norm' // &
      & '       f_norm        delu_norm     delux_norm    u_norm      cg_iter  attempt A|R' // &
      & '  bc_iter    solve', 'convergence header with the boundary-velocity iteration columns')
    ! Rows: (iteration, method/norm/CG index, solve).
    expect_rows = reshape([1, 1, 0, 2, 2, 0, 1, 1, 22], [3, 3])
    ok = .true.
    do i = 1_ik, 3_ik
      read (unit, *) ints, method, row(1:6), cgs(1), attempt_read, mark, bc_read
      j = int(expect_rows(2, i), ik)
      ok = ok .and. ints(1) == 5 .and. ints(2) == expect_rows(1, i)
      ok = ok .and. method == methods(j) .and. all(abs(row(1:6) - norms(:, j)) <= 0.0_rk)
      ok = ok .and. attempt_read == 1 .and. mark == 'A'
      ok = ok .and. bc_read(1) == 1 .and. bc_read(2) == expect_rows(3, i)
    end do
    read (unit, '(a)', iostat=ios) line
    call assert_true(ok .and. ios /= 0, 'convergence rows with the iteration and solve, exact')
    close (unit)
  end if

  ! Report, then end MPI.
  call test_summary(verbose=(my_rank == root_rank))
  call par_finalize()

contains

  !> Node test value from a global id and a step (fractions with 17 significant digits).
  !>
  !> Arguments:
  !>   id    (in)  integer(ik)  global node id
  !>   step  (in)  integer      step
  pure function node_value(id, step) result(v)

    integer(ik), intent(in) :: id
    integer, intent(in) :: step
    real(rk) :: v(3)

    v = [id/7.0_rk, -id/3.0_rk + step, 1.0e-200_rk*id]

  end function node_value

  !> Element test value (Mandel 6-vector) from a global id and a step.
  !>
  !> Arguments:
  !>   id    (in)  integer(ik)  global element id
  !>   step  (in)  integer      step
  !> Locals:
  !>   k     integer            component
  pure function elt_value(id, step) result(v)

    integer(ik), intent(in) :: id
    integer, intent(in) :: step
    real(rk) :: v(6)

    integer :: k

    v = [(1000.0_rk*k + id/3.0_rk + step, k = 1, 6)]

  end function elt_value

  !> Expected shear tensor components of elt_value (Mandel components divided by sqrt 2, as
  !> the writer computes them).
  !>
  !> Arguments:
  !>   id, step  (in)  global element id, step
  !> Locals:
  !>   v         real(rk)(6, 1)  the Mandel value
  pure function elt_value_shear(id, step) result(s)

    integer(ik), intent(in) :: id
    integer, intent(in) :: step
    real(rk) :: s(3)

    real(rk) :: v(6, 1)

    v(:, 1) = elt_value(id, step)
    v(4:6, :) = v(4:6, :)/sqrt2
    s = v(4:6, 1)

  end function elt_value_shear

  !> Number of blank-separated values on a line.
  !>
  !> Arguments:
  !>   text  (in)  character  the line
  !> Locals:
  !>   i     integer          character index
  !>   prev  logical          previous character is a blank
  pure function count_values(text) result(num)

    character(len=*), intent(in) :: text
    integer :: num

    integer :: i
    logical :: prev

    num = 0
    prev = .true.
    do i = 1, len_trim(text)
      if (prev .and. text(i:i) /= ' ') num = num + 1
      prev = text(i:i) == ' '
    end do

  end function count_values

  !> Whether two files have the same bytes.
  !>
  !> Arguments:
  !>   a, b  (in)  character  file names
  !> Locals:
  !>   ua, ub      integer          units
  !>   na, nb      integer(ik)      sizes
  !>   da, db      character(:)     contents
  function same_bytes(a, b) result(same)

    character(len=*), intent(in) :: a, b
    logical :: same

    integer :: ua, ub
    integer(ik) :: na, nb
    character(len=:), allocatable :: da, db

    open (newunit=ua, file=a, access='stream', form='unformatted', status='old', action='read')
    open (newunit=ub, file=b, access='stream', form='unformatted', status='old', action='read')
    inquire (unit=ua, size=na)
    inquire (unit=ub, size=nb)
    same = na == nb
    if (same) then
      allocate (character(len=na) :: da, db)
      read (ua) da
      read (ub) db
      same = da == db
    end if
    close (ua)
    close (ub)

  end function same_bytes

end program test_sim_writer
