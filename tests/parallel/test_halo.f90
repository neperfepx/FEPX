! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> MPI test of halo_mod: forward and reverse-add halo exchanges and distributed dot product and
!> norm; run on 1, 2, 3 and 4 ranks.
!>
!> Source : FEPX (ROADMAP step 2.4: "assembled nodal sums match the serial result exactly").
!>          Every rank parses the whole mesh as the serial reference. All values are integers or
!>          multiples of 1/4 small enough that every sum is exact in binary floating point, so the
!>          distributed results must equal the serial ones bit for bit, whatever the summation
!>          order.
!>
!> Arguments (command line): 1. fixture directory (tests/io/fixtures); 2. regression mesh
!> directory (tests/regression/cases/meshes).
!>
!> Checks, for poly_part4.msh (partition table, blocks) and tiny.msh (table, blocks; ranks without
!> elements on 3 and 4 ranks), with 1 and 3 values per node:
!>   - halo_forward: owned values v(g, c) = 4 g + c, ghosts -1 -> every local value is v; owned
!>     values unchanged; a second exchange right after (v + 1) is also exact;
!>   - halo_reverse_add: element contributions w(e, a, c) = e + a/4 + c/2 (element e, local node
!>     position a) accumulated on the local nodes -> owned values equal the serial sums over all
!>     elements; ghost values unchanged; a following halo_forward gives the sums on the ghosts;
!>   - dist_dot and dist_norm equal the serial values exactly.
!>
!> Locals:
!>   fixtures, meshes  character(512)     input directories
!>   ref               mesh_file_t        whole mesh (serial reference)
!>   mesh              mesh_t             this rank's mesh
!>   aux               character(512)(2)  auxiliary files of tiny.msh
!>   unit_log          integer            scratch log unit
!>   i_case            integer            case index
!>   num_comp          integer(ik)        values per node
program test_halo

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk, ik
  use log_mod, only: log_set_unit
  use parallel_mod, only: par_init, par_finalize, my_rank, root_rank
  use mesh_file_mod, only: mesh_file_t, mesh_file_load
  use mesh_mod, only: mesh_t
  use mesh_blocks_mod, only: mesh_load_blocks, mesh_load_table_or_blocks
  use halo_mod, only: halo_forward, halo_reverse_add, dist_dot, dist_norm
  use assert_mod, only: assert_true, assert_close, test_summary

  implicit none

  character(len=512) :: fixtures
  character(len=512) :: meshes
  type(mesh_file_t) :: ref
  type(mesh_t) :: mesh
  character(len=512) :: aux(2)
  integer :: unit_log
  integer :: i_case
  integer(ik) :: num_comp

  call par_init()
  call get_command_argument(1, fixtures)
  call get_command_argument(2, meshes)
  aux(1) = trim(fixtures) // '/tiny.opt'
  aux(2) = trim(fixtures) // '/tiny.ori'

  do i_case = 1, 4
    ! Load the case (reference and distributed mesh), with the log in a scratch file.
    open (newunit=unit_log, status='scratch', action='readwrite', form='formatted')
    call log_set_unit(unit_log)
    select case (i_case)
    case (1, 2)
      call mesh_file_load(ref, trim(meshes) // '/poly_part4.msh')
      if (i_case == 1) then
        call mesh_load_table_or_blocks(mesh, trim(meshes) // '/poly_part4.msh')
      else
        call mesh_load_blocks(mesh, trim(meshes) // '/poly_part4.msh')
      end if
    case (3, 4)
      call mesh_file_load(ref, trim(fixtures) // '/tiny.msh', aux)
      if (i_case == 3) then
        call mesh_load_table_or_blocks(mesh, trim(fixtures) // '/tiny.msh', aux)
      else
        call mesh_load_blocks(mesh, trim(fixtures) // '/tiny.msh', aux)
      end if
    end select
    call log_set_unit(output_unit)
    close (unit_log)

    ! Scalar and 3-component nodal vectors.
    do num_comp = 1_ik, 3_ik, 2_ik
      call check_halo(ref, mesh, num_comp, 'case ' // achar(48 + i_case) // ' comp ' // &
        & achar(48 + num_comp))
    end do
  end do

  call test_summary(verbose=(my_rank == root_rank))
  call par_finalize()

contains

  !> The checks of the program header for one mesh and one number of values per node.
  !>
  !> Arguments:
  !>   ref       (in)  type(mesh_file_t)  whole mesh
  !>   mesh      (in)  type(mesh_t)       this rank's mesh
  !>   num_comp  (in)  integer(ik)        values per node
  !>   label     (in)  character          prefix of the check names
  !> Locals:
  !>   x, y        real(rk)(:)     distributed nodal vectors
  !>   y_ghost     real(rk)(:)     ghost values of y before the reverse-add
  !>   serial      real(rk)(:)     serial assembled vector (all nodes)
  !>   expected    real(rk)        expected value
  !>   dot_serial  real(rk)        serial dot product
  !>   n_own       integer(ik)     owned entries (num_comp num_nodes_owned)
  !>   i, a, c, e, g  integer(ik)  local node, element node position, component, global
  !>                               element, global node
  !>   ok          logical         loop check result
  subroutine check_halo(ref, mesh, num_comp, label)

    type(mesh_file_t), intent(in) :: ref
    type(mesh_t), intent(in) :: mesh
    integer(ik), intent(in) :: num_comp
    character(len=*), intent(in) :: label

    real(rk), allocatable :: x(:)
    real(rk), allocatable :: y(:)
    real(rk), allocatable :: y_ghost(:)
    real(rk), allocatable :: serial(:)
    real(rk) :: expected
    real(rk) :: dot_serial
    integer(ik) :: n_own
    integer(ik) :: i
    integer(ik) :: a
    integer(ik) :: c
    integer(ik) :: e
    integer(ik) :: g
    logical :: ok

    n_own = num_comp*mesh%num_nodes_owned
    allocate (x(num_comp*mesh%num_nodes_loc), y(num_comp*mesh%num_nodes_loc))

    ! --- Forward: owned values v(g, c) = 4 g + c, ghosts -1 ----------------------------------

    x = -1.0_rk
    do i = 1_ik, mesh%num_nodes_owned
      do c = 1_ik, num_comp
        x(num_comp*(i - 1_ik) + c) = value_v(mesh%node_id_glob(i), c)
      end do
    end do
    call halo_forward(mesh, x, num_comp)
    call assert_true(all_local_equal(mesh, x, num_comp, 0.0_rk), label // ' forward')

    ! A second exchange right after the first (v + 1).
    x(1:n_own) = x(1:n_own) + 1.0_rk
    call halo_forward(mesh, x, num_comp)
    call assert_true(all_local_equal(mesh, x, num_comp, 1.0_rk), label // ' forward again')

    ! --- Reverse-add: element contributions --------------------------------------------------

    y = 0.0_rk
    do e = 1_ik, mesh%num_elts
      do a = 1_ik, mesh%num_elt_nodes
        i = mesh%elt_nodes_loc(a, e)
        do c = 1_ik, num_comp
          y(num_comp*(i - 1_ik) + c) = y(num_comp*(i - 1_ik) + c) + &
            & value_w(mesh%elt_id_glob(e), a, c)
        end do
      end do
    end do
    y_ghost = y(n_own + 1:)

    ! Serial assembly over all elements of the whole mesh.
    allocate (serial(num_comp*ref%num_nodes))
    serial = 0.0_rk
    do e = 1_ik, ref%num_elts
      do a = 1_ik, ref%num_elt_nodes
        g = ref%elt_nodes(a, e)
        do c = 1_ik, num_comp
          serial(num_comp*(g - 1_ik) + c) = serial(num_comp*(g - 1_ik) + c) + value_w(e, a, c)
        end do
      end do
    end do

    ! Owned values = serial sums; ghost values unchanged.
    call halo_reverse_add(mesh, y, num_comp)
    ok = .true.
    do i = 1_ik, mesh%num_nodes_owned
      g = mesh%node_id_glob(i)
      do c = 1_ik, num_comp
        if (abs(y(num_comp*(i - 1_ik) + c) - serial(num_comp*(g - 1_ik) + c)) > 0.0_rk) &
          & ok = .false.
      end do
    end do
    call assert_true(ok, label // ' reverse-add owned = serial')
    call assert_true(all(abs(y(n_own + 1:) - y_ghost) <= 0.0_rk), &
      & label // ' reverse-add ghosts unchanged')

    ! Forward: the ghosts get the assembled values.
    call halo_forward(mesh, y, num_comp)
    ok = .true.
    do i = 1_ik, mesh%num_nodes_loc
      g = mesh%node_id_glob(i)
      do c = 1_ik, num_comp
        if (abs(y(num_comp*(i - 1_ik) + c) - serial(num_comp*(g - 1_ik) + c)) > 0.0_rk) &
          & ok = .false.
      end do
    end do
    call assert_true(ok, label // ' assembled values on all local nodes')

    ! --- Dot product and norm ----------------------------------------------------------------

    ! x holds v + 1 on every local node (forward above).
    dot_serial = 0.0_rk
    do g = 1_ik, ref%num_nodes
      do c = 1_ik, num_comp
        expected = value_v(g, c) + 1.0_rk
        dot_serial = dot_serial + expected*serial(num_comp*(g - 1_ik) + c)
      end do
    end do
    call assert_close(dist_dot(x, y, n_own), dot_serial, label // ' dist_dot', abs_tol=0.0_rk, &
      & rel_tol=0.0_rk)
    dot_serial = 0.0_rk
    do g = 1_ik, ref%num_nodes
      do c = 1_ik, num_comp
        dot_serial = dot_serial + (value_v(g, c) + 1.0_rk)**2
      end do
    end do
    call assert_close(dist_norm(x, n_own), sqrt(dot_serial), label // ' dist_norm', &
      & abs_tol=0.0_rk, rel_tol=0.0_rk)

  end subroutine check_halo

  !> Whether every local value of a forwarded vector is v(g, c) + shift.
  !>
  !> Arguments:
  !>   mesh      (in)  type(mesh_t)  this rank's mesh
  !>   x         (in)  real(rk)(:)   nodal vector
  !>   num_comp  (in)  integer(ik)   values per node
  !>   shift     (in)  real(rk)      added to v
  !>   ok        (out, result) logical
  !> Locals:
  !>   i, c  integer(ik)  local node, component
  function all_local_equal(mesh, x, num_comp, shift) result(ok)

    type(mesh_t), intent(in) :: mesh
    real(rk), intent(in) :: x(:)
    integer(ik), intent(in) :: num_comp
    real(rk), intent(in) :: shift
    logical :: ok

    integer(ik) :: i
    integer(ik) :: c

    ok = .true.
    do i = 1_ik, mesh%num_nodes_loc
      do c = 1_ik, num_comp
        if (abs(x(num_comp*(i - 1_ik) + c) - (value_v(mesh%node_id_glob(i), c) + shift)) > &
          & 0.0_rk) ok = .false.
      end do
    end do

  end function all_local_equal

  !> Test value of component c of global node g: 4 g + c.
  !>
  !> Arguments:
  !>   g, c   (in)  integer(ik)  global node, component
  !>   v      (out, result) real(rk)
  pure function value_v(g, c) result(v)

    integer(ik), intent(in) :: g
    integer(ik), intent(in) :: c
    real(rk) :: v

    v = real(4_ik*g + c, rk)

  end function value_v

  !> Element contribution of global element e to its node at position a, component c:
  !> e + a/4 + c/2 (exact in binary).
  !>
  !> Arguments:
  !>   e, a, c  (in)  integer(ik)  global element, node position, component
  !>   w        (out, result) real(rk)
  pure function value_w(e, a, c) result(w)

    integer(ik), intent(in) :: e
    integer(ik), intent(in) :: a
    integer(ik), intent(in) :: c
    real(rk) :: w

    w = real(e, rk) + 0.25_rk*real(a, rk) + 0.5_rk*real(c, rk)

  end function value_w

end program test_halo
