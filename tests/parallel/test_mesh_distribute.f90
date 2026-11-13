! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> MPI test of mesh_distribute_mod (and par_scatterv): run on 1, 2, 3 and 4 ranks.
!>
!> Source : FEPX (ROADMAP step 2.2). Every rank also parses the mesh itself (mesh_file_mod) as
!>          the reference for its part.
!>
!> Arguments (command line): 1. fixture directory (tests/io/fixtures); 2. regression mesh
!> directory (tests/regression/cases/meshes).
!>
!> Checks:
!>   - even_block against v2's mre_decomp1d results (10 items on 4 ranks: 3 3 2 2; 2 items on 4
!>     ranks: 1 1 0 0) and its partition properties;
!>   - par_scatterv: parts of different sizes, including empty ones;
!>   - poly_part4.msh on 4 ranks: mesh_load uses the table (elements and owned nodes by
!>     partition); on 1 rank it takes the whole mesh; on 2 and 3 ranks, blocks (mesh_blocks_mod;
!>     mesh_load's error is tested by test_fepx_partition_*);
!>   - sx2_part4.msh on 4 ranks, table: the rank owning node 248 has no element containing it
!>     (Neper's layout; poly_part4.msh had such a node, 945, until its Neper 5.0.1-731-devel
!>     mesh, step 9.20);
!>   - a table that is not sequential (D108): the blocks reversed, by mesh_distribute;
!>   - poly_part4.msh in blocks (mesh_blocks_mod);
!>   - tiny.msh + tiny.opt + tiny.ori, table and blocks (ranks without elements on 3 and 4
!>     ranks): per-element values and element orientations of the local elements, per-elset
!>     values replicated;
!>   - for every distribution: each element on exactly one rank and each node owned by exactly
!>     one rank (sums over ranks); local lists increasing; ghosts = nodes of the local elements
!>     not owned, with their owners; coordinates, nodes, elsets, phases, rot_c2s of the local
!>     elements; nsets and fasets restricted to the local nodes and elements (local node numbers
!>     since step 2.3, mapped back to global ids here); periodicity and whole-mesh data
!>     replicated.
!>
!> Locals:
!>   fixtures    character(512)    fixture directory
!>   meshes      character(512)    regression mesh directory
!>   ref         mesh_file_t       reference (whole mesh parsed on this rank)
!>   mesh        mesh_t            this rank's distributed mesh
!>   aux         character(512)(2) auxiliary file names
!>   log_text    character(:)      captured log of this rank
!>   unit_log    integer           scratch log unit
!>   bounds      integer(ik)(2)    block bounds
!>   sizes       integer(ik)(4)    block sizes
!>   i_part      integer(ik)       block index
!>   send        integer(ik)(:)    scatter send buffer (root)
!>   counts      integer(ik)(:)    scatter counts (root)
!>   recv        integer(ik)(:)    received part
!>   send_r      real(rk)(:)       real scatter send buffer (root)
!>   recv_r      real(rk)(:)       received real part
!>   i           integer(ik)       index
!>   ok          logical           loop check result
program test_mesh_distribute

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk, ik
  use log_mod, only: log_set_unit
  use parallel_mod, only: par_init, par_finalize, par_scatterv, par_sum, num_ranks, my_rank, &
    & root_rank
  use mesh_file_mod, only: mesh_file_t, mesh_file_load
  use mesh_mod, only: mesh_t
  use mesh_distribute_mod, only: mesh_load, mesh_distribute
  use mesh_blocks_mod, only: even_block, mesh_load_blocks, mesh_load_table_or_blocks
  use assert_mod, only: assert_equal, assert_close, assert_true, test_summary

  implicit none

  character(len=512) :: fixtures
  character(len=512) :: meshes
  type(mesh_file_t) :: ref
  type(mesh_t) :: mesh
  character(len=512) :: aux(2)
  character(len=:), allocatable :: log_text
  integer :: unit_log
  integer(ik) :: bounds(2)
  integer(ik) :: sizes(4)
  integer(ik) :: i_part
  integer(ik), allocatable :: send(:)
  integer(ik), allocatable :: counts(:)
  integer(ik), allocatable :: recv(:)
  real(rk), allocatable :: send_r(:)
  real(rk), allocatable :: recv_r(:)
  integer(ik) :: i
  logical :: ok

  call par_init()
  call get_command_argument(1, fixtures)
  call get_command_argument(2, meshes)

  ! --- even_block (v2 mre_decomp1d) ------------------------------------------------------------

  do i_part = 0_ik, 3_ik
    bounds = even_block(10_ik, 4_ik, i_part)
    sizes(i_part + 1) = bounds(2) - bounds(1) + 1_ik
    if (i_part > 0_ik) call assert_equal(bounds(1), sum(sizes(1:i_part)) + 1_ik, &
      & 'even_block contiguous')
  end do
  call assert_true(all(sizes == [3_ik, 3_ik, 2_ik, 2_ik]), 'even_block 10 on 4')
  do i_part = 0_ik, 3_ik
    bounds = even_block(2_ik, 4_ik, i_part)
    sizes(i_part + 1) = bounds(2) - bounds(1) + 1_ik
  end do
  call assert_true(all(sizes == [1_ik, 1_ik, 0_ik, 0_ik]), 'even_block 2 on 4')

  ! --- par_scatterv: rank r receives r + 1 values (rank 0 of a 3+ rank run: also an empty part
  ! test via r = 1 receiving 0 below) ----------------------------------------------------------

  if (my_rank == root_rank) then
    counts = [(merge(0_ik, i + 1_ik, i == 1_ik), i = 0_ik, num_ranks - 1_ik)]
    send = [(i, i = 1_ik, sum(counts))]
    send_r = 0.5_rk*real(send, rk)
  else
    allocate (counts(0), send(0), send_r(0))
  end if
  call par_scatterv(send, counts, recv)
  call par_scatterv(send_r, counts, recv_r)
  call assert_equal(int(size(recv), ik), merge(0_ik, my_rank + 1_ik, my_rank == 1_ik), &
    & 'scatterv size')
  ok = .true.
  do i = 1_ik, int(size(recv), ik)
    ! Offset of this rank's part: sum of the earlier counts.
    if (recv(i) /= offset_before(my_rank) + i) ok = .false.
    if (abs(recv_r(i) - 0.5_rk*real(recv(i), rk)) > 0.0_rk) ok = .false.
  end do
  call assert_true(ok, 'scatterv values')

  ! --- poly_part4.msh, partition table (blocks on 2-3 ranks) -----------------------------------

  call mesh_file_load_quiet(ref, trim(meshes) // '/poly_part4.msh')
  call capture_start(unit_log)
  if (num_ranks == 1_ik .or. num_ranks == 4_ik) then
    call mesh_load(mesh, trim(meshes) // '/poly_part4.msh')
  else
    call mesh_load_table_or_blocks(mesh, trim(meshes) // '/poly_part4.msh')
  end if
  log_text = capture_end(unit_log)
  ! (The warning on 1 rank is fepx's setup summary's, from the mesh file's partition count.)
  if (num_ranks == 1_ik .or. num_ranks == 4_ik) &
    & call assert_equal(mesh%num_parts_file, 4_ik, 'poly4 partitions of the file')
  if (num_ranks == 4_ik) then
    call check_distribution(ref, mesh, 'poly4 table', .true.)
  else
    call check_distribution(ref, mesh, 'poly4 whole or blocks', .false.)
  end if

  ! --- sx2_part4.msh on 4 ranks, partition table: an owned node without a local element --------

  ! Neper's layout: node 248 is owned by rank 1 (partition 2) without a local element.
  if (num_ranks == 4_ik) then
    call mesh_file_load_quiet(ref, trim(meshes) // '/sx2_part4.msh')
    call capture_start(unit_log)
    call mesh_load(mesh, trim(meshes) // '/sx2_part4.msh')
    log_text = capture_end(unit_log)
    call check_distribution(ref, mesh, 'sx2_4 table', .true.)
    if (my_rank == 1_ik) then
      call assert_true(any(mesh%node_id_glob(1:mesh%num_nodes_owned) == 248_ik), &
        & 'sx2_4 rank 1 owns 248')
      call assert_true(.not. any(mesh%elt_nodes_glob == 248_ik), &
        & 'sx2_4 248 in no local element')
    end if
  end if

  ! --- poly_part4.msh, a table that is not sequential (D108): the blocks reversed ---------------

  call mesh_file_load_quiet(ref, trim(meshes) // '/poly_part4.msh')
  ref%num_parts = num_ranks
  do i_part = 0_ik, num_ranks - 1_ik
    bounds = even_block(ref%num_elts, num_ranks, i_part)
    ref%elt_part(bounds(1):bounds(2)) = i_part + 1_ik
    bounds = even_block(ref%num_nodes, num_ranks, i_part)
    ref%node_part(bounds(1):bounds(2)) = i_part + 1_ik
  end do
  ref%elt_part = ref%elt_part(ref%num_elts:1_ik:-1_ik)
  ref%node_part = ref%node_part(ref%num_nodes:1_ik:-1_ik)
  call capture_start(unit_log)
  call mesh_distribute(ref, mesh)
  log_text = capture_end(unit_log)
  call check_distribution(ref, mesh, 'poly4 not sequential', .true.)

  ! --- poly_part4.msh, blocks -------------------------------------------------------------------

  call mesh_file_load_quiet(ref, trim(meshes) // '/poly_part4.msh')
  call capture_start(unit_log)
  call mesh_load_blocks(mesh, trim(meshes) // '/poly_part4.msh')
  log_text = capture_end(unit_log)
  call check_distribution(ref, mesh, 'poly4 blocks', .false.)

  ! --- tiny.msh with auxiliary files, table and blocks ------------------------------------------

  aux(1) = trim(fixtures) // '/tiny.opt'
  aux(2) = trim(fixtures) // '/tiny.ori'
  call mesh_file_load_quiet(ref, trim(fixtures) // '/tiny.msh', aux)
  do i_part = 1_ik, 2_ik
    call capture_start(unit_log)
    if (i_part == 1_ik) then
      call mesh_load_table_or_blocks(mesh, trim(fixtures) // '/tiny.msh', aux)
    else
      call mesh_load_blocks(mesh, trim(fixtures) // '/tiny.msh', aux)
    end if
    log_text = capture_end(unit_log)
    call check_distribution(ref, mesh, 'tiny', i_part == 1_ik .and. num_ranks == 2_ik)
    ! Per-elset values replicated; per-element values of the local elements.
    call assert_true(mesh%values(1)%defined, 'tiny ElsetCrss defined')
    call assert_close(mesh%values(1)%values(:, 1), [100.0_rk, 110.0_rk], 'tiny ElsetCrss')
    call assert_true(mesh%values(5)%defined, 'tiny ElementCrssSat defined')
    ok = .true.
    do i = 1_ik, mesh%num_elts
      if (abs(mesh%values(5)%values(1, i) - ref%values(5)%values(1, mesh%elt_id_glob(i))) &
        & > 0.0_rk) ok = .false.
      if (mesh%values(5)%num_values(i) /= ref%values(5)%num_values(mesh%elt_id_glob(i))) &
        & ok = .false.
    end do
    call assert_true(ok, 'tiny ElementCrssSat local values')
    call assert_true(mesh%ori_from_elts, 'tiny element orientations')
  end do

  ! Report (summary on the root only) and end MPI.
  call test_summary(verbose=(my_rank == root_rank))
  call par_finalize()

contains

  !> Check a distributed mesh against the whole-mesh reference (see the program header).
  !>
  !> Arguments:
  !>   ref    (in)  type(mesh_file_t)  whole mesh parsed on this rank
  !>   mesh   (in)  type(mesh_t)       this rank's mesh
  !>   label  (in)  character          prefix of the check names
  !>   use_table  (in)  logical  the ranks' items are those of ref's table, not the blocks
  !> Locals:
  !>   elt_count   integer(ik)(:)  number of ranks holding each element (summed over ranks)
  !>   node_count  integer(ik)(:)  number of ranks owning each node (summed over ranks)
  !>   local_part  integer(ik)(:)  this rank's contribution
  !>   owner       integer(ik)(:)  owner rank of each node
  !>   bounds      integer(ik)(2)  block of this rank
  !>   ghosts      integer(ik)(:)  expected ghost nodes
  !>   mark        logical(:)      nodes of the local elements
  !>   i_set       integer         set index
  !>   expected    integer(ik)(:)  expected entries
  !>   is_local    logical(:)      local flag of each entry
  !>   e           integer(ik)     element number
  !>   i           integer(ik)     index
  !>   ok          logical         loop check result
  subroutine check_distribution(ref, mesh, label, use_table)

    type(mesh_file_t), intent(in) :: ref
    type(mesh_t), intent(in) :: mesh
    character(len=*), intent(in) :: label
    logical, intent(in) :: use_table

    integer(ik), allocatable :: elt_count(:)
    integer(ik), allocatable :: node_count(:)
    integer(ik), allocatable :: local_part(:)
    integer(ik), allocatable :: owner(:)
    integer(ik) :: bounds(2)
    integer(ik), allocatable :: ghosts(:)
    logical, allocatable :: mark(:)
    integer :: i_set
    integer(ik), allocatable :: expected(:)
    logical, allocatable :: is_local(:)
    integer(ik) :: e
    integer(ik) :: i
    logical :: ok

    ! Whole-mesh data.
    call assert_equal(mesh%num_nodes_glob, ref%num_nodes, label // ' num_nodes_glob')
    call assert_equal(mesh%num_elts_glob, ref%num_elts, label // ' num_elts_glob')
    call assert_equal(mesh%num_elt_nodes, ref%num_elt_nodes, label // ' num_elt_nodes')
    call assert_equal(mesh%num_elsets, ref%num_elsets, label // ' num_elsets')
    call assert_true(all(mesh%elset_id == ref%elset_id), label // ' elset_id')
    call assert_true(all(mesh%elset_phase == ref%elset_phase), label // ' elset_phase')
    call assert_equal(mesh%ori_descriptor, ref%ori_descriptor, label // ' ori_descriptor')
    call assert_equal(mesh%num_periodicity, ref%num_periodicity, label // ' num_periodicity')
    call assert_true(all(mesh%per_secondary == ref%per_secondary) .and. &
      & all(mesh%per_primary == ref%per_primary) .and. all(mesh%per_vect == ref%per_vect), &
      & label // ' periodicity')

    ! Each element on exactly one rank; each node owned by exactly one rank.
    allocate (local_part(ref%num_elts))
    local_part = 0_ik
    local_part(mesh%elt_id_glob) = 1_ik
    allocate (elt_count(ref%num_elts))
    call par_sum(local_part, elt_count)
    call assert_true(all(elt_count == 1_ik), label // ' each element on one rank')
    deallocate (local_part)
    allocate (local_part(ref%num_nodes), node_count(ref%num_nodes))
    local_part = 0_ik
    local_part(mesh%node_id_glob(1:mesh%num_nodes_owned)) = 1_ik
    call par_sum(local_part, node_count)
    call assert_true(all(node_count == 1_ik), label // ' each node owned once')

    ! Owner of every node (gathered as a sum of rank-tagged contributions).
    local_part = 0_ik
    local_part(mesh%node_id_glob(1:mesh%num_nodes_owned)) = my_rank
    allocate (owner(ref%num_nodes))
    call par_sum(local_part, owner)

    ! The rank's items: of the table, or of the blocks when num_ranks differs from its partitions
    ! (1 rank: the one block is the whole mesh).
    if (use_table) then
      call assert_true(all(ref%elt_part(mesh%elt_id_glob) == my_rank + 1_ik), &
        & label // ' elements of the partition')
      call assert_true(all(owner == ref%node_part - 1_ik), label // ' owners = node partitions')
    else
      bounds = even_block(ref%num_elts, num_ranks, my_rank)
      call assert_true(size(mesh%elt_id_glob) == max(0_ik, bounds(2) - bounds(1) + 1_ik), &
        & label // ' element block size')
      if (mesh%num_elts > 0_ik) call assert_true(mesh%elt_id_glob(1) == bounds(1), &
        & label // ' element block start')
      bounds = even_block(ref%num_nodes, num_ranks, my_rank)
      call assert_true(mesh%num_nodes_owned == max(0_ik, bounds(2) - bounds(1) + 1_ik), &
        & label // ' node block size')
    end if

    ! Local lists increasing.
    call assert_true(increasing(mesh%elt_id_glob), label // ' elements increasing')
    call assert_true(increasing(mesh%node_id_glob(1:mesh%num_nodes_owned)), &
      & label // ' owned increasing')
    call assert_true(increasing(mesh%node_id_glob(mesh%num_nodes_owned + 1:)), &
      & label // ' ghosts increasing')

    ! Ghosts: nodes of the local elements not owned here, with their owners.
    allocate (mark(ref%num_nodes))
    mark = .false.
    do i = 1_ik, mesh%num_elts
      mark(ref%elt_nodes(:, mesh%elt_id_glob(i))) = .true.
    end do
    ghosts = pack([(i, i = 1_ik, ref%num_nodes)], mark .and. owner /= my_rank)
    call assert_equal(mesh%num_nodes_ghost, int(size(ghosts), ik), label // ' ghost count')
    if (mesh%num_nodes_ghost == int(size(ghosts), ik)) then
      call assert_true(all(mesh%node_id_glob(mesh%num_nodes_owned + 1:) == ghosts), &
        & label // ' ghost ids')
      call assert_true(all(mesh%ghost_owner == owner(ghosts)), label // ' ghost owners')
    end if

    ! Coordinates of the local nodes; data of the local elements.
    call assert_true(all(abs(mesh%coo_s - ref%coo(:, mesh%node_id_glob)) <= 0.0_rk), &
      & label // ' coordinates')
    ok = .true.
    do i = 1_ik, mesh%num_elts
      e = mesh%elt_id_glob(i)
      if (any(mesh%elt_nodes_glob(:, i) /= ref%elt_nodes(:, e))) ok = .false.
      if (mesh%elt_elset(i) /= ref%elt_elset(e)) ok = .false.
      if (mesh%elt_phase(i) /= ref%elt_phase(e)) ok = .false.
      if (any(abs(mesh%rot_c2s(:, :, i) - ref%elt_rot_c2s(:, :, e)) > 0.0_rk)) ok = .false.
    end do
    call assert_true(ok, label // ' element data')

    ! Node sets restricted to the local nodes (file order).
    call assert_equal(size(mesh%nsets), size(ref%nsets), label // ' nset count')
    do i_set = 1, size(ref%nsets)
      is_local = [(any(mesh%node_id_glob == ref%nsets(i_set)%nodes(i)), &
        & i = 1_ik, int(size(ref%nsets(i_set)%nodes), ik))]
      expected = pack(ref%nsets(i_set)%nodes, is_local)
      call assert_true(mesh%nsets(i_set)%label == ref%nsets(i_set)%label .and. &
        & size(mesh%nsets(i_set)%nodes) == size(expected), label // ' nset size')
      if (size(mesh%nsets(i_set)%nodes) == size(expected)) then
        ! (local node numbers since step 2.3)
        call assert_true(all(mesh%node_id_glob(mesh%nsets(i_set)%nodes) == expected), &
          & label // ' nset nodes')
      end if
    end do

    ! Face sets restricted to the local elements (local element index).
    call assert_equal(size(mesh%fasets), size(ref%fasets), label // ' faset count')
    do i_set = 1, size(ref%fasets)
      is_local = [(any(mesh%elt_id_glob == ref%fasets(i_set)%elts(i)), &
        & i = 1_ik, int(size(ref%fasets(i_set)%elts), ik))]
      expected = pack(ref%fasets(i_set)%elts, is_local)
      call assert_true(size(mesh%fasets(i_set)%elts) == size(expected), label // ' faset size')
      if (size(mesh%fasets(i_set)%elts) == size(expected)) then
        call assert_true(all(mesh%elt_id_glob(mesh%fasets(i_set)%elts) == expected), &
          & label // ' faset elements')
        ! (local node numbers since step 2.3; compared flattened, in array element order)
        call assert_true(all(mesh%node_id_glob(reshape(mesh%fasets(i_set)%nodes, &
          & [size(mesh%fasets(i_set)%nodes)])) == reshape(ref%fasets(i_set)%nodes(:, &
          & pack([(i, i = 1_ik, int(size(ref%fasets(i_set)%elts), ik))], is_local)), &
          & [size(mesh%fasets(i_set)%nodes)])), label // ' faset nodes')
      end if
    end do

  end subroutine check_distribution

  !> Whether an integer array is strictly increasing.
  !>
  !> Arguments:
  !>   values  (in)  integer(ik)(:)  the array
  !>   inc     (out, result) logical
  pure function increasing(values) result(inc)

    integer(ik), intent(in) :: values(:)
    logical :: inc

    inc = .true.
    if (size(values) > 1) inc = all(values(2:) > values(:size(values) - 1))

  end function increasing

  !> Sum of the scatter counts of the ranks before a rank (counts r + 1, rank 1: 0).
  !>
  !> Arguments:
  !>   rank    (in)  integer(ik)  rank
  !>   offset  (out, result) integer(ik)
  !> Locals:
  !>   r  integer(ik)  earlier rank
  pure function offset_before(rank) result(offset)

    integer(ik), intent(in) :: rank
    integer(ik) :: offset

    integer(ik) :: r

    offset = 0_ik
    do r = 0_ik, rank - 1_ik
      if (r /= 1_ik) offset = offset + r + 1_ik
    end do

  end function offset_before

  !> Parse a mesh with the log redirected to a scratch file (reference on every rank).
  !>
  !> Arguments:
  !>   ref        (out)           type(mesh_file_t)  parsed mesh
  !>   name       (in)            character          mesh file name
  !>   aux_names  (in, optional)  character(*)(:)    auxiliary files
  !> Locals:
  !>   unit   integer  scratch unit (its content is discarded)
  subroutine mesh_file_load_quiet(ref, name, aux_names)

    type(mesh_file_t), intent(out) :: ref
    character(len=*), intent(in) :: name
    character(len=*), intent(in), optional :: aux_names(:)

    integer :: unit

    call capture_start(unit)
    call mesh_file_load(ref, name, aux_names)
    call log_set_unit(output_unit)
    close (unit)

  end subroutine mesh_file_load_quiet

  !> Start capturing this rank's log in a scratch file.
  !>
  !> Arguments:
  !>   unit  (out)  integer  scratch unit
  subroutine capture_start(unit)

    integer, intent(out) :: unit

    open (newunit=unit, status='scratch', action='readwrite', form='formatted')
    call log_set_unit(unit)

  end subroutine capture_start

  !> Stop capturing and return the log (lines joined with newlines); the root echoes it.
  !>
  !> Arguments:
  !>   unit  (in)  integer  scratch unit
  !>   text  (out, result) character(:), allocatable
  !> Locals:
  !>   line  character(1024)  one line
  !>   ios   integer          iostat of the reads
  function capture_end(unit) result(text)

    integer, intent(in) :: unit
    character(len=:), allocatable :: text

    character(len=1024) :: line
    integer :: ios

    call log_set_unit(output_unit)
    text = ''
    rewind (unit)
    do
      read (unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      text = text // trim(line) // new_line('a')
      if (my_rank == root_rank) write (output_unit, '(a)') trim(line)
    end do
    close (unit)

  end function capture_end

end program test_mesh_distribute
