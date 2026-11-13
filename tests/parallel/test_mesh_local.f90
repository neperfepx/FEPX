! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> MPI test of mesh_local_mod (local numbering, global <-> local maps, communication graph) and
!> par_alltoall / par_alltoallv: run on 1, 2, 3 and 4 ranks.
!>
!> Source : FEPX (ROADMAP step 2.3).
!>
!> Arguments (command line): 1. fixture directory (tests/io/fixtures); 2. regression mesh
!> directory (tests/regression/cases/meshes).
!>
!> Checks, for poly_part4.msh (partition table and blocks) and tiny.msh (ranks
!> without elements on 3 and 4 ranks):
!>   - maps: node_loc_from_glob(node_id_glob(i)) = i for every local node, and 0 for every
!>     global node that is not local; the same for the elements; node_dofs = 3(i - 1) + 1..3;
!>   - elt_nodes_loc maps back to elt_nodes_glob; nset and faset nodes are local numbers;
!>   - graph structure: neighbors increasing and not this rank; offsets increasing; the receive
!>     lists hold every ghost exactly once, each group owned by its neighbor, increasing global
!>     id within a group; the send lists hold owned nodes, increasing global id within a group;
!>   - graph consistency between ranks: each rank sends the global ids of its send lists to its
!>     neighbors (par_alltoallv); every rank receives exactly the global ids of its receive
!>     lists, group by group; the neighbor relation is symmetric (par_alltoall of flags).
!>
!> Locals:
!>   fixtures  character(512)     fixture directory
!>   meshes    character(512)     regression mesh directory
!>   mesh      mesh_t             this rank's mesh
!>   aux       character(512)(2)  auxiliary files of tiny.msh
!>   unit_log  integer            scratch log unit
!>   i_case    integer            case index
program test_mesh_local

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: ik
  use log_mod, only: log_set_unit
  use parallel_mod, only: par_init, par_finalize, par_alltoall, par_alltoallv, num_ranks, &
    & my_rank, root_rank
  use mesh_mod, only: mesh_t, node_loc_from_glob, &
    & elt_loc_from_glob, node_dofs
  use mesh_blocks_mod, only: mesh_load_blocks, mesh_load_table_or_blocks
  use assert_mod, only: assert_equal, assert_true, test_summary

  implicit none

  character(len=512) :: fixtures
  character(len=512) :: meshes
  type(mesh_t) :: mesh
  character(len=512) :: aux(2)
  integer :: unit_log
  integer :: i_case

  call par_init()
  call get_command_argument(1, fixtures)
  call get_command_argument(2, meshes)
  aux(1) = trim(fixtures) // '/tiny.opt'
  aux(2) = trim(fixtures) // '/tiny.ori'

  ! poly_part4 (table, even) and tiny (table, even); logs to a scratch file.
  do i_case = 1, 4
    open (newunit=unit_log, status='scratch', action='readwrite', form='formatted')
    call log_set_unit(unit_log)
    select case (i_case)
    case (1)
      call mesh_load_table_or_blocks(mesh, trim(meshes) // '/poly_part4.msh')
    case (2)
      call mesh_load_blocks(mesh, trim(meshes) // '/poly_part4.msh')
    case (3)
      call mesh_load_table_or_blocks(mesh, trim(fixtures) // '/tiny.msh', aux)
    case (4)
      call mesh_load_blocks(mesh, trim(fixtures) // '/tiny.msh', aux)
    end select
    call log_set_unit(output_unit)
    close (unit_log)
    call check_local(mesh, 'case ' // achar(48 + i_case))
  end do

  call test_summary(verbose=(my_rank == root_rank))
  call par_finalize()

contains

  !> Check the local numbering, the maps and the communication graph (program header).
  !>
  !> Arguments:
  !>   mesh   (in)  type(mesh_t)  this rank's mesh
  !>   label  (in)  character     prefix of the check names
  !> Locals:
  !>   is_local    logical(:)      whether each global node is local here
  !>   is_loc_elt  logical(:)      whether each global element is local here
  !>   g           integer(ik)     global id
  !>   i           integer(ik)     local index
  !>   a           integer(ik)     element node
  !>   k           integer(ik)     neighbor index
  !>   ok          logical         loop check result
  !>   covered     integer(ik)(:)  number of receive-list entries of each ghost
  !>   ids         integer(ik)(:)  global ids of a list
  !>   send_ids    integer(ik)(:)  global ids of the send lists, grouped by neighbor rank
  !>   send_cnt    integer(ik)(num_ranks)  number of ids sent to each rank
  !>   recv_cnt    integer(ik)(num_ranks)  number of ids received from each rank
  !>   recv_ids    integer(ik)(:)  ids received
  !>   expect_cnt  integer(ik)(num_ranks)  receive-list size of each rank
  !>   expect_ids  integer(ik)(:)  global ids of the receive lists, grouped by neighbor rank
  !>   flags       integer(ik)(num_ranks)  1 for each neighbor rank
  !>   flags_back  integer(ik)(num_ranks)  the neighbors' flags for this rank
  subroutine check_local(mesh, label)

    type(mesh_t), intent(in) :: mesh
    character(len=*), intent(in) :: label

    logical, allocatable :: is_local(:)
    logical, allocatable :: is_loc_elt(:)
    integer(ik) :: g
    integer(ik) :: i
    integer(ik) :: a
    integer(ik) :: k
    logical :: ok
    integer(ik), allocatable :: covered(:)
    integer(ik), allocatable :: ids(:)
    integer(ik), allocatable :: send_ids(:)
    integer(ik) :: send_cnt(num_ranks)
    integer(ik) :: recv_cnt(num_ranks)
    integer(ik), allocatable :: recv_ids(:)
    integer(ik) :: expect_cnt(num_ranks)
    integer(ik), allocatable :: expect_ids(:)
    integer(ik) :: flags(num_ranks)
    integer(ik) :: flags_back(num_ranks)

    ! --- Maps --------------------------------------------------------------------------------

    call assert_equal(mesh%num_nodes_loc, mesh%num_nodes_owned + mesh%num_nodes_ghost, &
      & label // ' num_nodes_loc')
    ok = .true.
    do i = 1_ik, mesh%num_nodes_loc
      if (node_loc_from_glob(mesh, mesh%node_id_glob(i)) /= i) ok = .false.
    end do
    call assert_true(ok, label // ' node map round trip')
    allocate (is_local(mesh%num_nodes_glob))
    is_local = .false.
    is_local(mesh%node_id_glob) = .true.
    ok = .true.
    do g = 1_ik, mesh%num_nodes_glob
      if ((node_loc_from_glob(mesh, g) > 0_ik) .neqv. is_local(g)) ok = .false.
    end do
    call assert_true(ok, label // ' node map of non-local nodes is 0')
    allocate (is_loc_elt(mesh%num_elts_glob))
    is_loc_elt = .false.
    is_loc_elt(mesh%elt_id_glob) = .true.
    ok = .true.
    do g = 1_ik, mesh%num_elts_glob
      i = elt_loc_from_glob(mesh, g)
      if ((i > 0_ik) .neqv. is_loc_elt(g)) ok = .false.
      if (i > 0_ik) then
        if (mesh%elt_id_glob(i) /= g) ok = .false.
      end if
    end do
    call assert_true(ok, label // ' element map')
    call assert_true(all(node_dofs(5_ik) == [13_ik, 14_ik, 15_ik]), label // ' node_dofs')

    ! --- Local numbers -----------------------------------------------------------------------

    call assert_true(all(shape(mesh%elt_nodes_loc) == shape(mesh%elt_nodes_glob)), &
      & label // ' elt_nodes_loc shape')
    ok = .true.
    do i = 1_ik, mesh%num_elts
      if (any(mesh%node_id_glob(mesh%elt_nodes_loc(:, i)) /= mesh%elt_nodes_glob(:, i))) &
        & ok = .false.
    end do
    call assert_true(ok, label // ' elt_nodes_loc maps to elt_nodes_glob')

    ! Element dof table: node a, component k at 3 (elt_nodes_loc(a) - 1) + k.
    ok = all(shape(mesh%elt_dofs_loc) == [3_ik*mesh%num_elt_nodes, mesh%num_elts])
    do i = 1_ik, mesh%num_elts
      do a = 1_ik, mesh%num_elt_nodes
        if (any(mesh%elt_dofs_loc(3*a - 2:3*a, i) /= 3_ik*(mesh%elt_nodes_loc(a, i) - 1_ik) + &
          & [1_ik, 2_ik, 3_ik])) ok = .false.
      end do
    end do
    call assert_true(ok, label // ' elt_dofs_loc from elt_nodes_loc')
    ok = .true.
    do k = 1_ik, int(size(mesh%nsets), ik)
      if (any(mesh%nsets(k)%nodes < 1_ik) .or. any(mesh%nsets(k)%nodes > mesh%num_nodes_loc)) &
        & ok = .false.
    end do
    do k = 1_ik, int(size(mesh%fasets), ik)
      if (any(mesh%fasets(k)%nodes < 1_ik) .or. &
        & any(mesh%fasets(k)%nodes > mesh%num_nodes_loc)) ok = .false.
    end do
    call assert_true(ok, label // ' set nodes are local numbers')

    ! --- Graph structure ---------------------------------------------------------------------

    associate (comm => mesh%comm)
      call assert_true(all(comm%neighbor /= my_rank) .and. increasing(comm%neighbor), &
        & label // ' neighbors increasing, not self')
      call assert_true(comm%recv_offset(1) == 1_ik .and. comm%send_offset(1) == 1_ik .and. &
        & all(comm%recv_offset(2:) >= comm%recv_offset(:comm%num_neighbors)) .and. &
        & all(comm%send_offset(2:) >= comm%send_offset(:comm%num_neighbors)), &
        & label // ' offsets')
      call assert_equal(comm%recv_offset(comm%num_neighbors + 1_ik) - 1_ik, &
        & mesh%num_nodes_ghost, label // ' receive list size')

      ! Every ghost received once, from its owner, in increasing global order per group.
      allocate (covered(mesh%num_nodes_loc))
      covered = 0_ik
      ok = .true.
      do k = 1_ik, comm%num_neighbors
        ids = comm%recv_nodes(comm%recv_offset(k):comm%recv_offset(k + 1_ik) - 1_ik)
        covered(ids) = covered(ids) + 1_ik
        if (any(ids <= mesh%num_nodes_owned)) ok = .false.
        if (any(mesh%ghost_owner(ids - mesh%num_nodes_owned) /= comm%neighbor(k))) ok = .false.
        if (.not. increasing(mesh%node_id_glob(ids))) ok = .false.
      end do
      call assert_true(ok, label // ' receive groups')
      call assert_true(all(covered(mesh%num_nodes_owned + 1:) == 1_ik), &
        & label // ' each ghost received once')

      ! Send lists: owned nodes, increasing global order per group.
      ok = .true.
      do k = 1_ik, comm%num_neighbors
        ids = comm%send_nodes(comm%send_offset(k):comm%send_offset(k + 1_ik) - 1_ik)
        if (any(ids < 1_ik) .or. any(ids > mesh%num_nodes_owned)) ok = .false.
        if (.not. increasing(mesh%node_id_glob(ids))) ok = .false.
      end do
      call assert_true(ok, label // ' send groups')

      ! --- Graph consistency between ranks ----------------------------------------------------

      ! Send the global ids of each send list to its neighbor.
      send_cnt = 0_ik
      expect_cnt = 0_ik
      allocate (send_ids(0), expect_ids(0))
      do k = 1_ik, comm%num_neighbors
        send_cnt(comm%neighbor(k) + 1_ik) = comm%send_offset(k + 1_ik) - comm%send_offset(k)
        expect_cnt(comm%neighbor(k) + 1_ik) = comm%recv_offset(k + 1_ik) - comm%recv_offset(k)
        send_ids = [send_ids, mesh%node_id_glob(comm%send_nodes(comm%send_offset(k): &
          & comm%send_offset(k + 1_ik) - 1_ik))]
        expect_ids = [expect_ids, mesh%node_id_glob(comm%recv_nodes(comm%recv_offset(k): &
          & comm%recv_offset(k + 1_ik) - 1_ik))]
      end do
      call par_alltoall(send_cnt, recv_cnt)
      call assert_true(all(recv_cnt == expect_cnt), label // ' counts sent = receive lists')
      call par_alltoallv(send_ids, send_cnt, recv_ids, recv_cnt)
      call assert_equal(int(size(recv_ids)), int(size(expect_ids)), label // ' ids size')
      if (size(recv_ids) == size(expect_ids)) then
        call assert_true(all(recv_ids == expect_ids), label // ' ids sent = receive lists')
      end if

      ! Symmetric neighbor relation.
      flags = 0_ik
      flags(comm%neighbor + 1_ik) = 1_ik
      call par_alltoall(flags, flags_back)
      call assert_true(all(flags == flags_back), label // ' neighbor relation symmetric')
    end associate

  end subroutine check_local

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

end program test_mesh_local
