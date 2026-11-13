! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> The .sim output directory (io_format.md §5; v2/Neper layout): node vectors (coo, disp, vel),
!> element variables (symmetric tensors such as stress and strain_el, orientations, slip-system
!> values, scalars and vectors; the list is elt_var_names), the face-set force tables (forces)
!> and the convergence record of the global iterations (convergence, step 6.5b).
!>
!> Source : port of v2 write_res_mod1.f90 (write_res, write_force_file_headers,
!>          write_force_file_data), write_res_mod2.f90 (write_res_vector, write_res_tensor,
!>          write_dot_sim_file_*) and units_mod.f90::open_output_files (directory layout, copies
!>          of the inputs). Contains no numerical method. Differences from v2 (D42, dev/
!>          cfg_keywords.md §7): every real is written with the one format sim_real_format
!>          (es24.16, an exact double round trip; v2: e13.7 after rounding to output_precision,
!>          forces e15.5); gathered values are ordered by global node or element id whatever the
!>          partitioning (v2 assumed contiguous blocks); an existing output directory is an
!>          error (v2 failed on its force files); no .sim summary file is written (v2
!>          write_dot_sim_file_*; Neper reads the directory without it, D109).
!>
!> Layout (dir, e.g. simulation.sim):
!>   dir/inputs/<file>                      copies of the input files
!>   dir/results/nodes/<var>/<var>.step<N>  one line per node (global id order), 3 values
!>   dir/results/elts/<var>/<var>.step<N>   one line per element (global id order), 6 tensor
!>                                          components 11 22 33 23 13 12, sample frame; ori: the
!>                                          parameters of the input's descriptor and convention
!>                                          (v2 write_res, CHANGES.md C3)
!>   dir/results/forces/<faset>             v2's header, then one row per call of
!>                                          sim_write_forces: step incr Fx Fy Fz area time
!>   dir/results/convergence                v2's post.conv header, then one row per global
!>                                          iteration (sim_write_convergence): incr iter
!>                                          SA|NR r_norm rx_norm f_norm delu_norm delux_norm
!>                                          u_norm cg_iter (step 6.5b; v2 writes post.conv in
!>                                          the working directory)
!>   dir/results/macro                      periodic conditions (step 8.3b, addendum §3.3,
!>                                          Q-P11): a header, then one row per call of
!>                                          sim_write_macro: step incr time, Fbar (11 12 13 21
!>                                          22 23 31 32 33), Fbar-dot (same order), and the
!>                                          volume-averaged Cauchy stress (11 22 33 23 13 12)
!> As in v2, forces and convergence are listed among the element results of the .sim file.
!>
!> Output formats (step 9.1, ROADMAP D93; the keyword output_format): ascii (the layout above,
!> the default), hdf5, or both. With hdf5 the root writes dir/simulation.h5 (serial HDF5,
!> hdf5_file_mod, D10) and its XDMF sidecar dir/simulation.xdmf; with hdf5 alone no text
!> result file is written (dir keeps inputs/, restart/, the two files).
!> The HDF5 layout (every dataset two-dimensional, (rows, columns) as the HDF5 tools show
!> it; the same names, values, orders, frames and units as the text files):
!>   /                                      attributes: format "fepx-hdf5 1", num_nodes,
!>                                          num_elts, num_elsets, num_ranks, element_type
!>                                          (tet10 or hex20), orientation ("<descriptor>:
!>                                          <convention>", or empty), inputs, node_results,
!>                                          elt_results (lists separated by blanks), num_steps
!>   /mesh/coo                (nodes, 3)    initial coordinates, sample frame [mm]
!>   /mesh/elts               (elts, 10|20) global node ids (1-based), FEPX's (Gmsh) order
!>   /mesh/elts_xdmf          (elts, 10|20) the same 0-based, in the XDMF (VTK) node order
!>   /mesh/elset, /mesh/phase (elts, 1)     elset id of the mesh file, phase
!>   /results/nodes/<var>/<var>.step<N>  (nodes, 3)
!>   /results/elts/<var>/<var>.step<N>   (elts, k); a slip-system variable has k = the
!>                                          largest number of systems, padded with NaN, and
!>                                          /results/elts/<var>/num_values (elts, 1) the
!>                                          number of each element
!>   /results/forces/<faset>  (rows, 7)     step incr Fx Fy Fz area time, one row per call
!>   /results/convergence     (rows, 14)    incr iter nr r_norm rx_norm f_norm delu_norm
!>                                          delux_norm u_norm cg_iter attempt accepted bc_iter
!>                                          solve (nr: 1 NR, 0 SA; accepted: 1 A, 0 R; the last
!>                                          four 1 1 0 0 where the text file has no column)
!>   /results/macro           (rows, 27)    the columns of results/macro
!>   /results/steps           (rows, 2)     step and time of every step written
!> Each table of rows has an attribute "columns" (the column names). The XDMF file (XDMF 3, a
!> temporal collection, one grid per step written, rewritten at every step) refers to the
!> datasets: topology Tet_10 or Hex_20 (/mesh/elts_xdmf), geometry the step's coo when printed
!> (else /mesh/coo), the node variables as vectors, the element variables as scalars (k = 1),
!> vectors (rotrate, rotrate_spin, rotrate_slip) or matrices of k columns (the text files'
!> order: symmetric tensors 11 22 33 23 13 12), and elset and phase.
!>
!> Parallel policy: every routine is collective; files are written by root_rank only, from data
!> gathered with par_gatherv. Errors on the root end every rank (parallel_mod abort handler).
!>
!> Public parameters:
!>   sim_real_format  edit descriptor of every real written (one constant, D42)
!> Public types:
!>   sim_writer_t     an open output directory (components at the type)
!> Public routines:
!>   sim_open          create the directory, copy the inputs, open the force files
!>   sim_write_nodes   write a node vector variable for one step
!>   sim_write_elts    write an element symmetric-tensor variable for one step
!>   sim_write_ori     write the element orientations for one step
!>   sim_write_elts_values  write an element variable of k values per element
!>   sim_write_elts_slip    write a slip-system element variable
!>   sim_write_forces  append one row to every force file
!>   sim_write_convergence  append the global iterations of an increment to the convergence
!>                     file
!>   sim_write_macro   append one row to the macroscopic file (periodic conditions)
!>   sim_close         write the HDF5 summary and close the output files
!>   sim_reopen        continue an existing directory from a restart point (step 8.4)
!>   sim_end_step      the end of a step's variables: its row of /results/steps and the XDMF
!>                     file (hdf5; step 9.1)
!>   sim_xdmf_node_order  the XDMF (VTK) order of an element's nodes (step 9.1)
!>   sim_restart_dir   the directory of the restart files, created when needed (step 8.4)
!>
!> Restart (step 8.4, dev/theory/addendum_restart.md §3 R3): the writer counts the lines of the
!> per-increment files (forces, convergence, macro; root); a restarted run reopens the
!> directory with sim_reopen, which removes the step files after the restart point, truncates
!> those files to their lines at the restart point and appends to them.
module sim_writer_mod

  use kinds_mod, only: rk, ik
  use constants_mod, only: sqrt2
  use error_mod, only: fatal_error, internal_error
  use parallel_mod, only: par_gatherv, par_barrier, num_ranks, my_rank, root_rank
  use mesh_mod, only: mesh_t
  use orientation_mod, only: ori_descriptor_name, ori_convention_name, ori_num_params, &
    & ori_max_params, descriptor_from_rot_c2s
  use hdf5_file_mod, only: h5_file_t, h5_create, h5_open, h5_close, h5_flush, h5_write_real, &
    & h5_write_int, h5_append_rows, h5_truncate_rows, h5_read_real, h5_dims, h5_exists, &
    & h5_delete, h5_write_attr_text, h5_write_attr_int
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

  implicit none

  private

  public :: sim_real_format, sim_writer_t
  public :: sim_open, sim_write_nodes, sim_write_elts, sim_write_elts_values, &
    & sim_write_elts_slip, sim_write_ori, sim_write_forces, sim_write_convergence, &
    & sim_write_macro, sim_close, sim_reopen, sim_restart_dir, sim_end_step, &
    & sim_xdmf_node_order

  !> Edit descriptor of every real in the result files (17 significant digits).
  character(len=*), parameter :: sim_real_format = 'es24.16'

  ! Node and element variables written so far, and v2's force-file header.
  character(len=*), parameter :: node_var_names = ' coo disp vel '
  character(len=*), parameter :: elt_var_names = ' stress strain_el ori crss slip sliprate ' // &
    & 'rss strain_pl strain_pl_eq defrate_pl spinrate rotrate work ' // &
    & 'strain strain_eq strain_el_eq stress_eq defrate defrate_eq defrate_pl_eq velgrad ' // &
    & 'rotrate_spin rotrate_slip work_pl workrate workrate_pl strain_pl_sample ' // &
    & 'strain_pl_sample_eq km_rho km_k2 km_rhorate dmg dmgrate '
  character(len=*), parameter :: force_header_1 = '% step     incr     Fx             ' // &
    & 'Fy             Fz             area a         time'
  character(len=*), parameter :: force_header_2 = '% -------- -------- -------------- ' // &
    & '-------------- -------------- -------------- --------------'
  ! HDF5 file, XDMF file, format attribute and the column names of the tables of rows (step 9.1).
  character(len=*), parameter :: h5_base = 'simulation.h5'
  character(len=*), parameter :: xdmf_base = 'simulation.xdmf'
  character(len=*), parameter :: h5_format = 'fepx-hdf5 1'
  character(len=*), parameter :: force_columns = 'step incr Fx Fy Fz area time'
  character(len=*), parameter :: conv_columns = 'incr iter nr r_norm rx_norm f_norm ' // &
    & 'delu_norm delux_norm u_norm cg_iter attempt accepted bc_iter solve'
  character(len=*), parameter :: macro_columns = 'step incr time Fbar11 Fbar12 Fbar13 ' // &
    & 'Fbar21 Fbar22 Fbar23 Fbar31 Fbar32 Fbar33 Fbardot11 Fbardot12 Fbardot13 Fbardot21 ' // &
    & 'Fbardot22 Fbardot23 Fbardot31 Fbardot32 Fbardot33 stress11 stress22 stress33 ' // &
    & 'stress23 stress13 stress12'
  ! Element variables written as XDMF vectors (sample-frame vectors).
  character(len=*), parameter :: xdmf_vector_names = ' rotrate rotrate_spin rotrate_slip '

  ! XDMF (VTK) node order of the elements (step 9.1): node i of the XDMF element is node
  ! order(i) of FEPX's (Gmsh's) order. VTK (VTK file formats, VTK_QUADRATIC_TETRA and
  ! VTK_QUADRATIC_HEXAHEDRON) puts the edge nodes of the tetrahedron on the edges (1,2) (2,3)
  ! (3,1) (1,4) (2,4) (3,4), those of the hexahedron on (1,2) (2,3) (3,4) (4,1) (5,6) (6,7)
  ! (7,8) (8,5) (1,5) (2,6) (3,7) (4,8) (1-based); FEPX's edge nodes are those of
  ! tet10_node_coo_ref and hex20_node_coo_ref (Gmsh's order). Checked by test_sim_writer_hdf5
  ! (each edge node at the middle of its edge).
  integer(ik), parameter :: tet10_xdmf_order(10) = [1, 2, 3, 4, 5, 6, 7, 8, 10, 9]
  integer(ik), parameter :: hex20_xdmf_order(20) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 14, 10, &
    & 17, 19, 20, 18, 11, 13, 15, 16]

  ! v2's convergence-file header (write_conv_file_headers).
  character(len=*), parameter :: conv_header = '%   incr     iter       NR     r_norm' // &
    & '        rx_norm       f_norm        delu_norm     delux_norm    u_norm      cg_iter'

  !> An open output directory.
  !>
  !> Components:
  !>   dir           character(:), allocatable  directory (e.g. "simulation.sim")
  !>   node_vars     character(32)(:)           node variables to write, print order
  !>   elt_vars      character(32)(:)           element variables (forces included), print order
  !>   inputs        character(256)(:)          base names of the copied inputs
  !>   faset_labels  character(64)(:)           face-set labels (force files), mesh order
  !>   force_units   integer(:)                 units of the open force files (root)
  !>   conv_unit     integer                    unit of the open convergence file (root; -1
  !>                                            when convergence is not printed)
  !>   conv_attempts logical                    the convergence file has the sub-incrementation
  !>                                            columns <attempt> <A|R> (Q35, D63)
  !>   conv_bc       logical                    the convergence file also has the columns
  !>                                            <bc_iter> <solve> of the boundary-velocity
  !>                                            iteration (step 8.2b)
  !>   has_macro     logical                    the macroscopic file is written (periodic
  !>                                            conditions)
  !>   macro_unit    integer                    unit of the open macroscopic file (root)
  !>   general       integer(ik)(4)             nodes, elements, elsets, ranks (HDF5 summary)
  !>   orides        character(:), allocatable  "<descriptor>:<convention>", or empty
  !>   num_steps     integer                    largest step written (HDF5 summary)
  !>   force_lines   integer                    lines of each force file, its header included
  !>                                            (root; the same for every face set; step 8.4)
  !>   conv_lines    integer                    lines of the convergence file (root)
  !>   macro_lines   integer                    lines of the macroscopic file (root)
  !>   reopened      logical                    the directory was continued (sim_reopen): the
  !>                                            .sim summary is replaced
  !>   ascii         logical                    the text files are written (output_format ascii)
  !>   hdf5          logical                    the HDF5 and XDMF files are written (output_format
  !>                                            hdf5; step 9.1)
  !>   h5            type(h5_file_t)            the open HDF5 file (root)
  !>   force_rows, conv_rows, macro_rows  integer  rows of the HDF5 tables of rows (root; each
  !>                                            face set the same)
  !>   h5_steps      integer(:)                 the steps written to the HDF5 file, in order (root)
  !>   h5_times      real(rk)(:)                their times [s] (root)
  !>   elt_type      integer(ik)                Gmsh type of the elements (11 or 17)
  !>   num_elt_nodes integer(ik)                nodes per element (10 or 20)
  type :: sim_writer_t
    character(len=:), allocatable :: dir
    character(len=32), allocatable :: node_vars(:)
    character(len=32), allocatable :: elt_vars(:)
    character(len=256), allocatable :: inputs(:)
    character(len=64), allocatable :: faset_labels(:)
    integer, allocatable :: force_units(:)
    integer :: conv_unit = -1
    logical :: conv_attempts = .false.
    logical :: conv_bc = .false.
    logical :: has_macro = .false.
    integer :: macro_unit = -1
    integer(ik) :: general(4) = 0_ik
    character(len=:), allocatable :: orides
    integer :: num_steps = 0
    integer :: force_lines = 0
    integer :: conv_lines = 0
    integer :: macro_lines = 0
    logical :: reopened = .false.
    logical :: ascii = .true.
    logical :: hdf5 = .false.
    type(h5_file_t) :: h5
    integer :: force_rows = 0
    integer :: conv_rows = 0
    integer :: macro_rows = 0
    integer, allocatable :: h5_steps(:)
    real(rk), allocatable :: h5_times(:)
    integer(ik) :: elt_type = 0_ik
    integer(ik) :: num_elt_nodes = 0_ik
  end type sim_writer_t

contains

  !> Create the output directory: the copies of the input files; with the text format the
  !> result directories of the printed variables and the force files with v2's header (when
  !> forces is printed); with the HDF5 format the HDF5 file with its summary attributes and the
  !> mesh (step 9.1). Collective.
  !>
  !> Arguments:
  !>   sim          (out) type(sim_writer_t)  the writer
  !>   dir          (in)  character           directory to create (must not exist)
  !>   print_vars   (in)  character(*)(:)     printed variables (config_t%print_vars)
  !>   mesh         (in)  type(mesh_t)        rank-local mesh (sizes, face sets, orientations)
  !>   input_files  (in)  character(*)(:)     input files to copy into dir/inputs (paths)
  !>   attempts     (in, optional) logical    the convergence file gets the columns <attempt>
  !>                                          <A|R> of sub-incrementation (Q35, D63; default
  !>                                          false: v2's columns)
  !>   bc_iterations (in, optional) logical   the convergence file also gets the columns
  !>                                          <bc_iter> <solve> of the boundary-velocity
  !>                                          iteration (step 8.2b; with attempts; default false)
  !>   macro        (in, optional) logical    write the macroscopic file results/macro
  !>                                          (periodic conditions, step 8.3b; default false)
  !>   ascii        (in, optional) logical    write the text files (default true)
  !>   hdf5         (in, optional) logical    write the HDF5 and XDMF files (default false)
  !> Locals:
  !>   i            integer                   index
  !>   exists       logical                   whether dir exists
  subroutine sim_open(sim, dir, print_vars, mesh, input_files, attempts, bc_iterations, &
    & macro, ascii, hdf5)

    type(sim_writer_t), intent(out) :: sim
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: print_vars(:)
    type(mesh_t), intent(in) :: mesh
    character(len=*), intent(in) :: input_files(:)
    logical, intent(in), optional :: macro
    logical, intent(in), optional :: attempts
    logical, intent(in), optional :: bc_iterations
    logical, intent(in), optional :: ascii, hdf5

    integer :: i
    logical :: exists

    ! The variables, the formats and the summary data (every rank).
    call init_writer(sim, dir, print_vars, mesh, input_files, macro, ascii, hdf5)
    if (present(attempts)) sim%conv_attempts = attempts
    if (present(bc_iterations)) sim%conv_bc = bc_iterations
    if (sim%conv_bc) sim%conv_attempts = .true.

    ! Directories (an existing output directory is an error: results are never mixed); the
    ! copies of the inputs (root).
    if (my_rank == root_rank) then
      inquire (file=dir // '/.', exist=exists)
      if (exists) call fatal_error('sim_open', 'output directory ' // dir // ' already exists')
      call make_dir(dir // '/inputs')
      do i = 1, size(input_files)
        call copy_file(trim(input_files(i)), dir // '/inputs/' // trim(sim%inputs(i)))
      end do
    end if

    ! HDF5 file: the summary, then the mesh (gathered from every rank; collective).
    if (sim%hdf5) then
      if (my_rank == root_rank) then
        call h5_create(sim%h5, dir // '/' // h5_base)
        call h5_write_summary(sim)
      end if
      call h5_write_mesh(sim, mesh)
    end if
    if (my_rank /= root_rank .or. .not. sim%ascii) return

    ! Text format: the result directories.
    call make_dir(dir // '/results/nodes')
    call make_dir(dir // '/results/elts')
    do i = 1, size(sim%node_vars)
      call make_dir(dir // '/results/nodes/' // trim(sim%node_vars(i)))
    end do
    do i = 1, size(sim%elt_vars)
      if (sim%elt_vars(i) /= 'forces' .and. sim%elt_vars(i) /= 'convergence') &
        & call make_dir(dir // '/results/elts/' // trim(sim%elt_vars(i)))
    end do

    ! Force files with v2's header.
    if (any(sim%elt_vars == 'forces')) then
      call make_dir(dir // '/results/forces')
      deallocate (sim%force_units)
      allocate (sim%force_units(size(sim%faset_labels)))
      do i = 1, size(sim%faset_labels)
        sim%force_units(i) = open_new(dir // '/results/forces/' // trim(sim%faset_labels(i)))
        write (sim%force_units(i), '(a)') force_header_1
        write (sim%force_units(i), '(a)') force_header_2
      end do
      sim%force_lines = 2
    end if

    ! Convergence file with v2's header (with sub-incrementation: two more columns; with the
    ! boundary-velocity iteration two more).
    if (any(sim%elt_vars == 'convergence')) then
      sim%conv_unit = open_new(dir // '/results/convergence')
      call write_conv_header(sim)
    end if

    ! Macroscopic file (periodic conditions, step 8.3b).
    if (sim%has_macro) then
      sim%macro_unit = open_new(dir // '/results/macro')
      call write_macro_header(sim)
    end if

  end subroutine sim_open

  !> Write the convergence file's header line (v2's, with the columns of sub-incrementation and
  !> of the boundary-velocity iteration when present); counts it. Root only.
  !>
  !> Arguments:
  !>   sim  (inout) type(sim_writer_t)  the writer (convergence file open)
  subroutine write_conv_header(sim)

    type(sim_writer_t), intent(inout) :: sim

    if (sim%conv_bc) then
      write (sim%conv_unit, '(a)') conv_header // '  attempt A|R  bc_iter    solve'
    else if (sim%conv_attempts) then
      write (sim%conv_unit, '(a)') conv_header // '  attempt A|R'
    else
      write (sim%conv_unit, '(a)') conv_header
    end if
    sim%conv_lines = sim%conv_lines + 1

  end subroutine write_conv_header

  !> Write the macroscopic file's header line; counts it. Root only.
  !>
  !> Arguments:
  !>   sim  (inout) type(sim_writer_t)  the writer (macroscopic file open)
  subroutine write_macro_header(sim)

    type(sim_writer_t), intent(inout) :: sim

    write (sim%macro_unit, '(a)') '% step incr time | Fbar 11 12 13 21 22 23 31 32 33 | ' // &
      & 'Fbar-dot (same order) [1/s] | stress (volume average) 11 22 33 23 13 12 [MPa]'
    sim%macro_lines = sim%macro_lines + 1

  end subroutine write_macro_header

  !> The writer's variables by entity, its formats and its summary data (every rank; sim_open,
  !> sim_reopen).
  !>
  !> Arguments:
  !>   sim, dir, print_vars, mesh, input_files, macro, ascii, hdf5  as sim_open
  !> Locals:
  !>   i    integer                    index
  !>   var  character(:), allocatable  a variable name
  subroutine init_writer(sim, dir, print_vars, mesh, input_files, macro, ascii, hdf5)

    type(sim_writer_t), intent(inout) :: sim
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: print_vars(:)
    type(mesh_t), intent(in) :: mesh
    character(len=*), intent(in) :: input_files(:)
    logical, intent(in), optional :: macro
    logical, intent(in), optional :: ascii, hdf5

    integer :: i
    character(len=:), allocatable :: var

    ! Variables by entity (forces with the element results, as v2).
    sim%dir = dir
    allocate (sim%node_vars(0), sim%elt_vars(0))
    do i = 1, size(print_vars)
      var = trim(print_vars(i))
      if (index(node_var_names, ' ' // var // ' ') > 0) then
        sim%node_vars = [character(len=32) :: sim%node_vars, var]
      else if (index(elt_var_names, ' ' // var // ' ') > 0 .or. var == 'forces' .or. &
        & var == 'convergence') then
        sim%elt_vars = [character(len=32) :: sim%elt_vars, var]
      else
        call internal_error('sim_open', 'output variable ' // var // ' is not written yet')
      end if
    end do

    ! Summary data: mesh sizes, ranks, orientation descriptor, inputs, face sets.
    sim%general = [mesh%num_nodes_glob, mesh%num_elts_glob, mesh%num_elsets, num_ranks]
    sim%orides = ''
    if (mesh%has_ori) sim%orides = ori_descriptor_name(mesh%ori_descriptor) // ':' // &
      & ori_convention_name(mesh%ori_convention)
    allocate (sim%inputs(size(input_files)))
    do i = 1, size(input_files)
      sim%inputs(i) = base_name(trim(input_files(i)))
    end do
    allocate (sim%faset_labels(size(mesh%fasets)), sim%force_units(0))
    do i = 1, size(mesh%fasets)
      sim%faset_labels(i) = mesh%fasets(i)%label
    end do
    if (present(macro)) sim%has_macro = macro

    ! Formats (at least one) and the element type (step 9.1).
    if (present(ascii)) sim%ascii = ascii
    if (present(hdf5)) sim%hdf5 = hdf5
    if (.not. (sim%ascii .or. sim%hdf5)) call internal_error('sim_open', 'no output format')
    sim%elt_type = mesh%elt_type
    sim%num_elt_nodes = mesh%num_elt_nodes
    allocate (sim%h5_steps(0), sim%h5_times(0))

  end subroutine init_writer

  !> Continue an existing output directory from a restart point (step 8.4, addendum_restart.md
  !> §3 R3): the inputs replaced by this run's; with the text format its directories completed
  !> for the printed variables, the step files of the step after the restart point removed
  !> (written by an interrupted run), the force, convergence and macroscopic files truncated to
  !> their lines at the restart point and appended to (a missing one created with its header;
  !> the convergence file gets a new header line when its columns change); with the HDF5 format
  !> (step 9.1) the same in the HDF5 file: the datasets of that step removed, the tables of
  !> rows cut back to their rows at the restart point, the summary and the XDMF file
  !> rewritten. Collective.
  !>
  !> Arguments:
  !>   sim, dir, print_vars, mesh, input_files, attempts, bc_iterations, macro  as sim_open
  !>   last_step    (in)  integer   largest step written at the restart point
  !>   force_lines, conv_lines, macro_lines  (in)  integer  lines of the force files, the
  !>                                convergence file and the macroscopic file at the restart
  !>                                point (0: the file was not written)
  !>   prev_attempts, prev_bc  (in)  logical  the convergence file's columns at the restart point
  !>   ascii, hdf5  (in)  logical   the formats (as sim_open; the restart point's: the keyword
  !>                                output_format is strict on a restart)
  !>   rows         (in)  integer(3) rows of the HDF5 force tables, convergence table and
  !>                                macroscopic table at the restart point
  !> Locals:
  !>   i        integer                    index
  !>   exists   logical                    whether dir (or the HDF5 file) exists
  !>   names    character(:), allocatable  every variable name the writer knows
  !>   var      character(:), allocatable  one of them
  !>   pos      integer                    position in names
  !>   step_txt character(16)              the step after the restart point, as text
  !>   steps    real(rk)(:, :)             (2, rows) the HDF5 table of the steps written
  !>   keep     integer                    its rows at the restart point
  subroutine sim_reopen(sim, dir, print_vars, mesh, input_files, attempts, bc_iterations, &
    & macro, last_step, force_lines, conv_lines, macro_lines, prev_attempts, prev_bc, ascii, &
    & hdf5, rows)

    type(sim_writer_t), intent(out) :: sim
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: print_vars(:)
    type(mesh_t), intent(in) :: mesh
    character(len=*), intent(in) :: input_files(:)
    logical, intent(in) :: attempts, bc_iterations, macro
    integer, intent(in) :: last_step, force_lines, conv_lines, macro_lines
    logical, intent(in) :: prev_attempts, prev_bc
    logical, intent(in) :: ascii, hdf5
    integer, intent(in) :: rows(3)

    integer :: i, pos, keep
    logical :: exists
    character(len=:), allocatable :: names, var
    character(len=16) :: step_txt
    real(rk), allocatable :: steps(:, :)

    ! The variables, the formats and the summary data (every rank); the steps written so far.
    call init_writer(sim, dir, print_vars, mesh, input_files, macro, ascii, hdf5)
    sim%reopened = .true.
    sim%num_steps = last_step
    sim%conv_attempts = attempts
    sim%conv_bc = bc_iterations
    if (sim%conv_bc) sim%conv_attempts = .true.
    if (my_rank /= root_rank) return

    ! The directory must exist; this run's inputs (they replace the copies of the previous
    ! run's).
    inquire (file=dir // '/.', exist=exists)
    if (.not. exists) call fatal_error('sim_reopen', 'restart on: output directory ' // dir // &
      & ' not found')
    call make_dir(dir // '/inputs')
    do i = 1, size(input_files)
      call copy_file(trim(input_files(i)), dir // '/inputs/' // trim(sim%inputs(i)), .true.)
    end do
    write (step_txt, '(i0)') last_step + 1

    ! HDF5 file (step 9.1): the step after the restart point removed, the tables of rows cut
    ! back, the steps written so far, the summary and the XDMF file.
    if (sim%hdf5) then
      inquire (file=dir // '/' // h5_base, exist=exists)
      if (.not. exists) call fatal_error('sim_reopen', 'restart on: ' // dir // '/' // &
        & h5_base // ' not found')
      call h5_open(sim%h5, dir // '/' // h5_base)
      call for_each_name(node_var_names, 'nodes')
      call for_each_name(elt_var_names, 'elts')
      sim%force_rows = rows(1)
      sim%conv_rows = rows(2)
      sim%macro_rows = rows(3)
      do i = 1, size(sim%faset_labels)
        if (h5_exists(sim%h5, '/results/forces/' // trim(sim%faset_labels(i)))) call &
          & h5_truncate_rows(sim%h5, '/results/forces/' // trim(sim%faset_labels(i)), rows(1))
      end do
      if (h5_exists(sim%h5, '/results/convergence')) call h5_truncate_rows(sim%h5, &
        & '/results/convergence', rows(2))
      if (h5_exists(sim%h5, '/results/macro')) call h5_truncate_rows(sim%h5, '/results/macro', &
        & rows(3))
      if (h5_exists(sim%h5, '/results/steps')) then
        call h5_read_real(sim%h5, '/results/steps', steps)
        keep = count(nint(steps(1, :)) <= last_step)
        call h5_truncate_rows(sim%h5, '/results/steps', keep)
        sim%h5_steps = nint(steps(1, 1:keep))
        sim%h5_times = steps(2, 1:keep)
      end if
      call h5_write_summary(sim)
      call write_xdmf(sim)
      call h5_flush(sim%h5)
    end if
    if (.not. sim%ascii) return

    ! Text format: the result directories of this run's variables.
    do i = 1, size(sim%node_vars)
      call make_dir(dir // '/results/nodes/' // trim(sim%node_vars(i)))
    end do
    do i = 1, size(sim%elt_vars)
      if (sim%elt_vars(i) /= 'forces' .and. sim%elt_vars(i) /= 'convergence') &
        & call make_dir(dir // '/results/elts/' // trim(sim%elt_vars(i)))
    end do

    ! The step files after the restart point (an interrupted run writes at most the step after
    ! it: every completed step writes a restart point), of every variable.
    names = node_var_names
    do while (len_trim(names) > 0)
      names = adjustl(names)
      pos = index(names, ' ')
      var = names(1:pos - 1)
      names = names(pos:)
      call delete_if_exists(dir // '/results/nodes/' // var // '/' // var // '.step' // &
        & trim(step_txt))
    end do
    names = elt_var_names
    do while (len_trim(names) > 0)
      names = adjustl(names)
      pos = index(names, ' ')
      var = names(1:pos - 1)
      names = names(pos:)
      call delete_if_exists(dir // '/results/elts/' // var // '/' // var // '.step' // &
        & trim(step_txt))
    end do

    ! Force files: truncated to the restart point, or created with v2's header.
    if (any(sim%elt_vars == 'forces')) then
      call make_dir(dir // '/results/forces')
      deallocate (sim%force_units)
      allocate (sim%force_units(size(sim%faset_labels)))
      do i = 1, size(sim%faset_labels)
        if (force_lines > 0) then
          sim%force_units(i) = open_truncated(dir // '/results/forces/' // &
            & trim(sim%faset_labels(i)), force_lines)
        else
          sim%force_units(i) = open_new(dir // '/results/forces/' // trim(sim%faset_labels(i)))
          write (sim%force_units(i), '(a)') force_header_1
          write (sim%force_units(i), '(a)') force_header_2
        end if
      end do
      sim%force_lines = max(force_lines, 2)
    end if

    ! Convergence file: truncated, a new header line when the columns change, or created.
    if (any(sim%elt_vars == 'convergence')) then
      if (conv_lines > 0) then
        sim%conv_unit = open_truncated(dir // '/results/convergence', conv_lines)
        sim%conv_lines = conv_lines
        if (sim%conv_attempts .neqv. prev_attempts .or. (sim%conv_bc .neqv. prev_bc)) &
          & call write_conv_header(sim)
      else
        sim%conv_unit = open_new(dir // '/results/convergence')
        call write_conv_header(sim)
      end if
    end if

    ! Macroscopic file: truncated, or created with its header.
    if (sim%has_macro) then
      if (macro_lines > 0) then
        sim%macro_unit = open_truncated(dir // '/results/macro', macro_lines)
        sim%macro_lines = macro_lines
      else
        sim%macro_unit = open_new(dir // '/results/macro')
        call write_macro_header(sim)
      end if
    end if

  contains

    !> Remove the HDF5 datasets of the step after the restart point of every variable of a
    !> list (names separated by blanks) of one entity.
    !>
    !> Arguments:
    !>   list    (in) character  the names
    !>   entity  (in) character  nodes or elts
    !> Locals:
    !>   rest    character(:), allocatable  the names not yet done
    !>   word    character(:), allocatable  one name
    !>   at      integer                    position of the blank after it
    subroutine for_each_name(list, entity)

      character(len=*), intent(in) :: list, entity

      character(len=:), allocatable :: rest, word
      integer :: at

      rest = list
      do while (len_trim(rest) > 0)
        rest = adjustl(rest)
        at = index(rest, ' ')
        word = rest(1:at - 1)
        rest = rest(at:)
        call h5_delete(sim%h5, '/results/' // entity // '/' // word // '/' // word // '.step' // &
          & trim(step_txt))
      end do

    end subroutine for_each_name

  end subroutine sim_reopen

  !> The directory of the restart files, <dir>/restart (step 8.4, addendum_restart.md §5),
  !> created by the root when needed. Collective (a barrier after the creation).
  !>
  !> Arguments:
  !>   sim   (in)  type(sim_writer_t)  the writer
  !>   path  (out, result) character(:), allocatable  the directory
  function sim_restart_dir(sim) result(path)

    type(sim_writer_t), intent(in) :: sim
    character(len=:), allocatable :: path

    path = sim%dir // '/restart'
    if (my_rank == root_rank) call make_dir(path)
    call par_barrier()

  end function sim_restart_dir

  !> Append one row to the macroscopic file (periodic conditions, step 8.3b, addendum §3.3,
  !> Q-P11): the step, the increment, the time, Fbar, Fbar-dot and the volume-averaged Cauchy
  !> stress. The values must be the global ones; the root writes them. Does nothing without the
  !> file. Collective in effect (no communication).
  !>
  !> Arguments:
  !>   sim     (in)  type(sim_writer_t)  the writer
  !>   step    (in)  integer             step (0: the initial state)
  !>   incr    (in)  integer             increment (cumulative)
  !>   time    (in)  real(rk)            time [s]
  !>   fbar    (in)  real(rk)(3, 3)      macroscopic deformation gradient [-]
  !>   fdot    (in)  real(rk)(3, 3)      its rate [1/s]
  !>   stress  (in)  real(rk)(6)         volume-averaged Cauchy stress, tensor components 11 22
  !>                                     33 23 13 12 [MPa]
  !> Locals:
  !>   row     real(rk)(27)              the row of the HDF5 table (the text file's columns)
  subroutine sim_write_macro(sim, step, incr, time, fbar, fdot, stress)

    type(sim_writer_t), intent(inout) :: sim
    integer, intent(in) :: step, incr
    real(rk), intent(in) :: time, fbar(3, 3), fdot(3, 3), stress(6)

    real(rk) :: row(27)

    if (my_rank /= root_rank .or. .not. sim%has_macro) return
    ! Row by row (11 12 13 21 ...): the transposes of the column-major arrays.
    if (sim%ascii) then
      write (sim%macro_unit, '(i8,1x,i8,25(1x,' // sim_real_format // '))') step, incr, time, &
        & transpose(fbar), transpose(fdot), stress
      sim%macro_lines = sim%macro_lines + 1
      flush (sim%macro_unit)
    end if
    if (sim%hdf5) then
      row = [real(step, rk), real(incr, rk), time, reshape(transpose(fbar), [9]), &
        & reshape(transpose(fdot), [9]), stress]
      call append_row(sim, '/results/macro', row, macro_columns, sim%macro_rows == 0)
      sim%macro_rows = sim%macro_rows + 1
    end if

  end subroutine sim_write_macro

  !> Write a node vector variable (e.g. coo, disp) for one step: the owned nodes of every rank,
  !> gathered and ordered by global id, one line of 3 values per node. Collective. Does nothing
  !> when the variable is not printed.
  !>
  !> Arguments:
  !>   sim     (inout) type(sim_writer_t)  the writer
  !>   mesh    (in)    type(mesh_t)        rank-local mesh
  !>   name    (in)    character           variable name
  !>   values  (in)    real(rk)(:, :)      (3, num_nodes_loc or num_nodes_owned) values of the
  !>                                       local nodes (owned first), sample frame
  !>   step    (in)    integer             step number (0: initial state)
  !> Locals:
  !>   n_own   integer(ik)                 owned nodes
  !>   ids     integer(ik)(:)              gathered global ids (root)
  !>   vals    real(rk)(:)                 gathered values, 3 per node (root)
  !>   table   real(rk)(:, :)              (3, num_nodes_glob) values by global id (root)
  !>   i       integer(ik)                 index
  !>   unit    integer                     file unit
  subroutine sim_write_nodes(sim, mesh, name, values, step)

    type(sim_writer_t), intent(inout) :: sim
    type(mesh_t), intent(in) :: mesh
    character(len=*), intent(in) :: name
    real(rk), intent(in) :: values(:, :)
    integer, intent(in) :: step

    integer(ik) :: n_own, i
    integer(ik), allocatable :: ids(:)
    real(rk), allocatable :: vals(:), table(:, :)
    integer :: unit

    if (.not. any(sim%node_vars == name)) return
    if (size(values, 1) /= 3 .or. size(values, 2) < mesh%num_nodes_owned) call internal_error( &
      & 'sim_write_nodes', name // ': expected 3 values per owned node')

    ! Owned nodes of every rank on the root, then placed by global id.
    n_own = mesh%num_nodes_owned
    call par_gatherv(mesh%node_id_glob(1:n_own), ids)
    call par_gatherv(reshape(values(:, 1:n_own), [3*n_own]), vals)
    sim%num_steps = max(sim%num_steps, step)
    if (my_rank /= root_rank) return
    allocate (table(3, mesh%num_nodes_glob))
    call place(ids, vals, table, 'sim_write_nodes', name)

    ! One line per node; one row per node of the HDF5 dataset.
    if (sim%ascii) then
      unit = open_new(sim%dir // '/results/nodes/' // name // '/' // name // '.step' // &
        & int_text(step))
      do i = 1_ik, mesh%num_nodes_glob
        write (unit, '(*(' // sim_real_format // '))') table(:, i)
      end do
      close (unit)
    end if
    if (sim%hdf5) then
      call h5_write_real(sim%h5, '/results/nodes/' // name // '/' // name // '.step' // &
        & int_text(step), table)
      call h5_flush(sim%h5)
    end if

  end subroutine sim_write_nodes

  !> Write an element symmetric-tensor variable (e.g. stress, strain_el) for one step: the
  !> elements of every rank, gathered and ordered by global id, one line of the 6 tensor
  !> components 11 22 33 23 13 12 per element (Mandel shear components divided by sqrt 2).
  !> Collective. Does nothing when the variable is not printed.
  !>
  !> Theory : Mandel 6-vector {a} = (a11, a22, a33, sqrt2 a23, sqrt2 a13, sqrt2 a12)
  !>          (conventions.md; recoupled_formulation.md §2).
  !>
  !> Arguments:
  !>   sim     (inout) type(sim_writer_t)  the writer
  !>   mesh    (in)    type(mesh_t)        rank-local mesh
  !>   name    (in)    character           variable name
  !>   values  (in)    real(rk)(:, :)      (6, num_elts) Mandel values of the local elements
  !>                                       (at the centroid point), sample frame
  !>   step    (in)    integer             step number (0: initial state)
  !> Locals:
  !>   comps   real(rk)(:, :)              (6, num_elts) tensor components
  subroutine sim_write_elts(sim, mesh, name, values, step)

    type(sim_writer_t), intent(inout) :: sim
    type(mesh_t), intent(in) :: mesh
    character(len=*), intent(in) :: name
    real(rk), intent(in) :: values(:, :)
    integer, intent(in) :: step

    real(rk), allocatable :: comps(:, :)

    if (.not. any(sim%elt_vars == name)) return
    if (size(values, 1) /= 6 .or. size(values, 2) /= mesh%num_elts) call internal_error( &
      & 'sim_write_elts', name // ': expected 6 values per local element')

    ! Tensor components: shear Mandel components divided by sqrt 2.
    comps = values
    comps(4:6, :) = comps(4:6, :)/sqrt2

    ! Gathered, ordered by global id and written.
    call write_elt_table(sim, mesh, name, comps, step, 'sim_write_elts')

  end subroutine sim_write_elts

  !> Write an element variable of k values per element as given (e.g. a scalar, k = 1, or a
  !> 3-vector already in its output order) for one step: the elements of every rank, gathered
  !> and ordered by global id, one line of k values per element. Collective. Does nothing when
  !> the variable is not printed.
  !>
  !> Arguments:
  !>   sim     (inout) type(sim_writer_t)  the writer
  !>   mesh    (in)    type(mesh_t)        rank-local mesh
  !>   name    (in)    character           variable name
  !>   values  (in)    real(rk)(:, :)      (k, num_elts) values of the local elements (at the
  !>                                       centroid point)
  !>   step    (in)    integer             step number (0: initial state)
  subroutine sim_write_elts_values(sim, mesh, name, values, step)

    type(sim_writer_t), intent(inout) :: sim
    type(mesh_t), intent(in) :: mesh
    character(len=*), intent(in) :: name
    real(rk), intent(in) :: values(:, :)
    integer, intent(in) :: step

    if (.not. any(sim%elt_vars == name)) return
    if (size(values, 2) /= mesh%num_elts) call internal_error('sim_write_elts_values', &
      & name // ': expected the values of every local element')

    ! Gathered, ordered by global id and written.
    call write_elt_table(sim, mesh, name, values, step, 'sim_write_elts_values')

  end subroutine sim_write_elts_values

  !> Write a slip-system variable (e.g. crss, sliprate) for one step: one line per element with
  !> the values of its phase's slip systems, in their order (the first num_values(i) rows of
  !> column i; v2 write_res writes each element's own number of systems). Collective. Does
  !> nothing when the variable is not printed.
  !>
  !> Arguments:
  !>   sim         (inout) type(sim_writer_t)  the writer
  !>   mesh        (in)    type(mesh_t)        rank-local mesh
  !>   name        (in)    character           variable name
  !>   values      (in)    real(rk)(:, :)      (num_slip_max, num_elts) values of the local
  !>                                           elements (at the centroid point)
  !>   num_values  (in)    integer(ik)(:)      (num_elts) slip systems of each local element
  !>   step        (in)    integer             step number (0: initial state)
  subroutine sim_write_elts_slip(sim, mesh, name, values, num_values, step)

    type(sim_writer_t), intent(inout) :: sim
    type(mesh_t), intent(in) :: mesh
    character(len=*), intent(in) :: name
    real(rk), intent(in) :: values(:, :)
    integer(ik), intent(in) :: num_values(:)
    integer, intent(in) :: step

    if (.not. any(sim%elt_vars == name)) return
    if (size(values, 2) /= mesh%num_elts .or. size(num_values) /= mesh%num_elts) &
      & call internal_error('sim_write_elts_slip', name // ': expected the values of every ' // &
      & 'local element')

    ! Gathered, ordered by global id and written with each element's count.
    call write_elt_table(sim, mesh, name, values, step, 'sim_write_elts_slip', num_values)

  end subroutine sim_write_elts_slip

  !> Write the element orientations for one step: the parameters of the input's descriptor in
  !> the input's convention (conventions §2, D11; v2 write_res, CHANGES.md C3), one line per
  !> element. Collective. Does nothing when ori is not printed.
  !>
  !> Arguments:
  !>   sim      (inout) type(sim_writer_t)  the writer
  !>   mesh     (in)    type(mesh_t)        rank-local mesh (descriptor and convention)
  !>   rot_c2s  (in)    real(rk)(:, :, :)   (3, 3, num_elts) crystal-to-sample rotations of the
  !>                                        local elements (at the centroid point)
  !>   step     (in)    integer             step number (0: initial state)
  !> Locals:
  !>   num      integer                     parameters of the descriptor
  !>   params   real(rk)(:, :)              (num, num_elts) parameters
  !>   all_params  real(rk)(ori_max_params)   parameters of one element (unused entries 0)
  !>   i        integer(ik)                 element
  subroutine sim_write_ori(sim, mesh, rot_c2s, step)

    type(sim_writer_t), intent(inout) :: sim
    type(mesh_t), intent(in) :: mesh
    real(rk), intent(in) :: rot_c2s(:, :, :)
    integer, intent(in) :: step

    integer :: num
    real(rk), allocatable :: params(:, :)
    real(rk) :: all_params(ori_max_params)
    integer(ik) :: i

    if (.not. any(sim%elt_vars == 'ori')) return
    if (.not. mesh%has_ori) call internal_error('sim_write_ori', 'the mesh has no orientations')
    if (size(rot_c2s, 1) /= 3 .or. size(rot_c2s, 2) /= 3 .or. size(rot_c2s, 3) /= &
      & mesh%num_elts) call internal_error('sim_write_ori', 'expected a 3 x 3 rotation per ' // &
      & 'local element')

    ! Descriptor parameters of each element.
    num = ori_num_params(mesh%ori_descriptor)
    allocate (params(num, mesh%num_elts))
    do i = 1_ik, mesh%num_elts
      all_params = descriptor_from_rot_c2s(mesh%ori_descriptor, mesh%ori_convention, &
        & rot_c2s(:, :, i))
      params(:, i) = all_params(1:num)
    end do

    ! Gathered, ordered by global id and written.
    call write_elt_table(sim, mesh, 'ori', params, step, 'sim_write_ori')

  end subroutine sim_write_ori

  !> Gather an element table (k values per local element) on the root, order it by global id
  !> and write it, one line per element, to results/elts/<name>/<name>.step<N>. Collective.
  !>
  !> Arguments:
  !>   sim      (inout) type(sim_writer_t)  the writer
  !>   mesh     (in)    type(mesh_t)        rank-local mesh
  !>   name     (in)    character           variable name
  !>   values   (in)    real(rk)(:, :)      (k, num_elts) values of the local elements
  !>   step     (in)    integer             step number
  !>   routine  (in)    character           calling routine, for messages
  !>   num_values (in, optional) integer(ik)(:)  (num_elts) values to write per element (the
  !>                                       first ones of its column; default all k)
  !> Locals:
  !>   k        integer                     values per element
  !>   ids      integer(ik)(:)              gathered global ids (root)
  !>   vals     real(rk)(:)                 gathered values, k per element (root)
  !>   table    real(rk)(:, :)              (k, num_elts_glob) values by global id (root)
  !>   counts   real(rk)(:)                 gathered num_values, as reals (root)
  !>   count_table  real(rk)(:, :)          (1, num_elts_glob) num_values by global id, or k
  !>                                        (root)
  !>   i        integer(ik)                 index
  !>   unit     integer                     file unit
  !>   path     character(:), allocatable   the HDF5 group of the variable
  subroutine write_elt_table(sim, mesh, name, values, step, routine, num_values)

    type(sim_writer_t), intent(inout) :: sim
    type(mesh_t), intent(in) :: mesh
    character(len=*), intent(in) :: name
    real(rk), intent(in) :: values(:, :)
    integer, intent(in) :: step
    character(len=*), intent(in) :: routine
    integer(ik), intent(in), optional :: num_values(:)

    integer :: k, unit
    integer(ik), allocatable :: ids(:)
    real(rk), allocatable :: vals(:), table(:, :), counts(:), count_table(:, :)
    integer(ik) :: i
    character(len=:), allocatable :: path

    ! Elements of every rank on the root, then placed by global id.
    k = size(values, 1)
    call par_gatherv(mesh%elt_id_glob, ids)
    call par_gatherv(reshape(values, [k*mesh%num_elts]), vals)
    if (present(num_values)) call par_gatherv(real(num_values, rk), counts)
    sim%num_steps = max(sim%num_steps, step)
    if (my_rank /= root_rank) return
    allocate (table(k, mesh%num_elts_glob))
    call place(ids, vals, table, routine, name)
    ! The values per element: num_values by global id, or k (always allocated).
    allocate (count_table(1, mesh%num_elts_glob))
    count_table = real(k, rk)
    if (present(num_values)) call place(ids, counts, count_table, routine, name)

    ! One line per element (its first num_values values when given).
    if (sim%ascii) then
      unit = open_new(sim%dir // '/results/elts/' // name // '/' // name // '.step' // &
        & int_text(step))
      do i = 1_ik, mesh%num_elts_glob
        if (present(num_values)) then
          write (unit, '(*(' // sim_real_format // '))') table(1:nint(count_table(1, i)), i)
        else
          write (unit, '(*(' // sim_real_format // '))') table(:, i)
        end if
      end do
      close (unit)
    end if

    ! One row per element of the HDF5 dataset; with num_values the values after each
    ! element's count are NaN, and the counts are written once (they do not change).
    if (sim%hdf5) then
      path = '/results/elts/' // name
      if (present(num_values)) then
        do i = 1_ik, mesh%num_elts_glob
          table(nint(count_table(1, i)) + 1:, i) = ieee_value(1.0_rk, ieee_quiet_nan)
        end do
        if (.not. h5_exists(sim%h5, path // '/num_values')) call h5_write_int(sim%h5, &
          & path // '/num_values', nint(count_table, ik))
      end if
      call h5_write_real(sim%h5, path // '/' // name // '.step' // int_text(step), table)
      call h5_flush(sim%h5)
    end if

  end subroutine write_elt_table

  !> Append one row to every force file: step, increment, the force on the face set, its area
  !> and the time (v2 write_force_file_data, uniaxial form). The values must be the global
  !> values (reduced over the ranks); the root writes them. Does nothing when forces is not
  !> printed.
  !>
  !> Arguments:
  !>   sim     (in)  type(sim_writer_t)  the writer
  !>   step    (in)  integer             step number
  !>   incr    (in)  integer             increment number (cumulative)
  !>   loads   (in)  real(rk)(:, :)      (3, num_fasets) force on each face set, sample frame [N]
  !>   areas   (in)  real(rk)(:)         (num_fasets) area of each face set [mm^2]
  !>   time    (in)  real(rk)            time [s]
  !> Locals:
  !>   i       integer                   face-set index
  !>   row     real(rk)(7)               the row of an HDF5 table (the text file's columns)
  subroutine sim_write_forces(sim, step, incr, loads, areas, time)

    type(sim_writer_t), intent(inout) :: sim
    integer, intent(in) :: step, incr
    real(rk), intent(in) :: loads(:, :)
    real(rk), intent(in) :: areas(:)
    real(rk), intent(in) :: time

    integer :: i
    real(rk) :: row(7)

    if (.not. any(sim%elt_vars == 'forces')) return
    if (size(loads, 1) /= 3 .or. size(loads, 2) /= size(sim%faset_labels) .or. &
      & size(areas) /= size(sim%faset_labels)) call internal_error('sim_write_forces', &
      & 'expected a force and an area per face set')
    if (my_rank /= root_rank) return

    ! One row per face set, in its text file and its HDF5 table.
    do i = 1, size(sim%faset_labels)
      if (sim%ascii) then
        write (sim%force_units(i), '(2i8,*(' // sim_real_format // '))') step, incr, &
          & loads(:, i), areas(i), time
        ! Flushed: a restart point counts the lines written (step 8.4).
        flush (sim%force_units(i))
      end if
      if (sim%hdf5) then
        row = [real(step, rk), real(incr, rk), loads(:, i), areas(i), time]
        call append_row(sim, '/results/forces/' // trim(sim%faset_labels(i)), row, &
          & force_columns, sim%force_rows == 0)
      end if
    end do
    if (sim%ascii) sim%force_lines = sim%force_lines + 1
    if (sim%hdf5) sim%force_rows = sim%force_rows + 1

  end subroutine sim_write_forces

  !> Append the global iterations of one increment to the convergence file: one row per
  !> iteration with the increment and iteration numbers, the method, v2's six norms and the CG
  !> iterations (v2 write_conv_file_data; format: sim_real_format for the norms). The values must
  !> be the global values (reduced over the ranks); the root writes them. Does nothing when
  !> convergence is not printed. With sub-incrementation (sim_open attempts, Q35, D63,
  !> addendum_subincrementation.md §5) every attempt of an increment is written under the
  !> increment's number, with two more columns: the attempt (1, 2, ... in the order solved) and
  !> A (accepted) or R (rejected). With the boundary-velocity iteration (sim_open bc_iterations,
  !> step 8.2b) every global solve of the iteration is written, with two more columns: the
  !> boundary-velocity iteration (1, 2, ...) and the solve: 0 at the iteration's velocities,
  !> 11, 22 or 33 for the solve with that component's face velocity perturbed.
  !>
  !> Arguments:
  !>   sim      (in)  type(sim_writer_t)  the writer
  !>   incr     (in)  integer             increment number (cumulative; 0: the initialization
  !>                                      increment)
  !>   methods  (in)  character(2)(:)     (num_iters) 'SA' or 'NR' of each iteration
  !>   norms    (in)  real(rk)(:, :)      (6, num_iters) r_norm, rx_norm, f_norm, delu_norm,
  !>                                      delux_norm, u_norm (increment_info_t%iter_norms)
  !>   cg_iters (in)  integer(:)          (num_iters) CG iterations of each iteration
  !>   attempt  (in, optional) integer    attempt within the increment (written only with
  !>                                      sim_open attempts; default 1)
  !>   accepted (in, optional) logical    the attempt was accepted (likewise; default true)
  !>   bc_iter  (in, optional) integer    boundary-velocity iteration (written only with sim_open
  !>                                      bc_iterations; default 0)
  !>   bc_solve (in, optional) integer    0, or the perturbed component 11, 22, 33 (likewise;
  !>                                      default 0)
  !> Locals:
  !>   i        integer                   iteration
  !>   num      integer                   the attempt
  !>   mark     character(1)              A or R
  !>   iter_bc, solve_bc  integer         the boundary-velocity iteration and solve
  !>   row      real(rk)(14)              a row of the HDF5 table
  subroutine sim_write_convergence(sim, incr, methods, norms, cg_iters, attempt, accepted, &
    & bc_iter, bc_solve)

    type(sim_writer_t), intent(inout) :: sim
    integer, intent(in) :: incr
    character(len=2), intent(in) :: methods(:)
    real(rk), intent(in) :: norms(:, :)
    integer, intent(in) :: cg_iters(:)
    integer, intent(in), optional :: attempt
    logical, intent(in), optional :: accepted
    integer, intent(in), optional :: bc_iter, bc_solve

    integer :: i, num, iter_bc, solve_bc
    character(len=1) :: mark
    real(rk) :: row(14)

    if (.not. any(sim%elt_vars == 'convergence')) return
    if (size(norms, 1) /= 6 .or. size(norms, 2) /= size(methods) .or. &
      & size(cg_iters) /= size(methods)) call internal_error('sim_write_convergence', &
      & 'expected six norms and a CG count per iteration')
    if (my_rank /= root_rank) return

    ! The attempt and its mark (sub-incrementation).
    num = 1
    if (present(attempt)) num = attempt
    mark = 'A'
    if (present(accepted)) then
      if (.not. accepted) mark = 'R'
    end if

    ! The boundary-velocity iteration and solve (step 8.2b).
    iter_bc = 0
    solve_bc = 0
    if (present(bc_iter)) iter_bc = bc_iter
    if (present(bc_solve)) solve_bc = bc_solve

    ! HDF5 table: one row per global iteration, every column (nr: 1 NR, 0 SA; accepted: 1, 0).
    if (sim%hdf5) then
      do i = 1, size(methods)
        row = [real(incr, rk), real(i, rk), merge(1.0_rk, 0.0_rk, methods(i) == 'NR'), &
          & norms(:, i), real(cg_iters(i), rk), real(num, rk), merge(1.0_rk, 0.0_rk, &
          & mark == 'A'), real(iter_bc, rk), real(solve_bc, rk)]
        call append_row(sim, '/results/convergence', row, conv_columns, sim%conv_rows == 0)
        sim%conv_rows = sim%conv_rows + 1
      end do
    end if
    if (.not. sim%ascii) return

    ! Text file: one row per global iteration.
    do i = 1, size(methods)
      if (sim%conv_bc) then
        write (sim%conv_unit, '(i8,1x,i8,1x,a2,6(1x,' // sim_real_format // &
          & '),1x,i8,1x,i8,1x,a1,1x,i8,1x,i8)') incr, i, methods(i), norms(:, i), cg_iters(i), &
          & num, mark, iter_bc, solve_bc
      else if (sim%conv_attempts) then
        write (sim%conv_unit, '(i8,1x,i8,1x,a2,6(1x,' // sim_real_format // &
          & '),1x,i8,1x,i8,1x,a1)') incr, i, methods(i), norms(:, i), cg_iters(i), num, mark
      else
        write (sim%conv_unit, '(i8,1x,i8,1x,a2,6(1x,' // sim_real_format // '),1x,i8)') incr, &
          & i, methods(i), norms(:, i), cg_iters(i)
      end if
    end do
    sim%conv_lines = sim%conv_lines + size(methods)
    ! Flushed: a failed increment ends the run before sim_close, and its rows show why.
    flush (sim%conv_unit)

  end subroutine sim_write_convergence

  !> Write the HDF5 summary and the XDMF file and close the HDF5 file (HDF5 format); close the
  !> force and convergence files (text format). Collective.
  !>
  !> Arguments:
  !>   sim        (inout) type(sim_writer_t)  the writer
  !>   num_steps  (in, optional) integer      completed steps (default: the largest step
  !>                                          written)
  !> Locals:
  !>   i          integer                     index
  subroutine sim_close(sim, num_steps)

    type(sim_writer_t), intent(inout) :: sim
    integer, intent(in), optional :: num_steps

    integer :: i

    if (present(num_steps)) sim%num_steps = num_steps
    if (my_rank /= root_rank) return

    ! HDF5 file (step 9.1): the summary (its steps), the XDMF file, closed.
    if (sim%hdf5) then
      call h5_write_summary(sim)
      call write_xdmf(sim)
      call h5_close(sim%h5)
    end if
    if (.not. sim%ascii) return

    ! Force files and the convergence file.
    do i = 1, size(sim%force_units)
      close (sim%force_units(i))
    end do
    if (sim%conv_unit >= 0) close (sim%conv_unit)
    if (sim%has_macro) close (sim%macro_unit)

  end subroutine sim_close

  !> The end of a step's variables (step 9.1): with the HDF5 format, the step's row of
  !> /results/steps (step and time) and the XDMF file rewritten with the step (so that the
  !> results of a run that stops can be opened). Nothing with the text format. Collective in
  !> effect (no communication; the root writes).
  !>
  !> Arguments:
  !>   sim   (inout) type(sim_writer_t)  the writer
  !>   step  (in)    integer             step number (0: initial state)
  !>   time  (in)    real(rk)            its time [s]
  subroutine sim_end_step(sim, step, time)

    type(sim_writer_t), intent(inout) :: sim
    integer, intent(in) :: step
    real(rk), intent(in) :: time

    if (my_rank /= root_rank .or. .not. sim%hdf5) return
    call append_row(sim, '/results/steps', [real(step, rk), time], 'step time', &
      & size(sim%h5_steps) == 0)
    sim%h5_steps = [sim%h5_steps, step]
    sim%h5_times = [sim%h5_times, time]
    call write_xdmf(sim)

  end subroutine sim_end_step

  !> Append one row to an HDF5 table of rows, its column names written with the first row;
  !> the file flushed. Root only.
  !>
  !> Arguments:
  !>   sim      (in) type(sim_writer_t)  the writer (HDF5 file open)
  !>   path     (in) character           the table
  !>   row      (in) real(rk)(:)         the row
  !>   columns  (in) character           the column names, separated by blanks
  !>   first    (in) logical             the table has no row yet
  subroutine append_row(sim, path, row, columns, first)

    type(sim_writer_t), intent(in) :: sim
    character(len=*), intent(in) :: path
    real(rk), intent(in) :: row(:)
    character(len=*), intent(in) :: columns
    logical, intent(in) :: first

    call h5_append_rows(sim%h5, path, reshape(row, [size(row), 1]))
    if (first) call h5_write_attr_text(sim%h5, path, 'columns', columns)
    call h5_flush(sim%h5)

  end subroutine append_row

  !> Write (or rewrite) the summary attributes of the HDF5 file's root group (step 9.1): the
  !> format, the sizes, the element type, the orientation descriptor, the inputs, the result
  !> lists and the completed steps. Root only.
  !>
  !> Arguments:
  !>   sim  (in) type(sim_writer_t)  the writer (HDF5 file open)
  !> Locals:
  !>   elt_name  character(5)        tet10 or hex20
  subroutine h5_write_summary(sim)

    type(sim_writer_t), intent(in) :: sim

    character(len=5) :: elt_name

    elt_name = 'tet10'
    if (sim%elt_type == 17_ik) elt_name = 'hex20'
    call h5_write_attr_text(sim%h5, '/', 'format', h5_format)
    call h5_write_attr_int(sim%h5, '/', 'num_nodes', [sim%general(1)])
    call h5_write_attr_int(sim%h5, '/', 'num_elts', [sim%general(2)])
    call h5_write_attr_int(sim%h5, '/', 'num_elsets', [sim%general(3)])
    call h5_write_attr_int(sim%h5, '/', 'num_ranks', [sim%general(4)])
    call h5_write_attr_text(sim%h5, '/', 'element_type', elt_name)
    call h5_write_attr_text(sim%h5, '/', 'orientation', sim%orides)
    call h5_write_attr_text(sim%h5, '/', 'inputs', list_text(sim%inputs))
    call h5_write_attr_text(sim%h5, '/', 'node_results', list_text(sim%node_vars))
    call h5_write_attr_text(sim%h5, '/', 'elt_results', list_text(sim%elt_vars))
    call h5_write_attr_int(sim%h5, '/', 'num_steps', [int(sim%num_steps, ik)])
    call h5_flush(sim%h5)

  end subroutine h5_write_summary

  !> Write the mesh to the HDF5 file (step 9.1): the initial coordinates of the nodes, the
  !> element nodes (global ids, FEPX's order; and 0-based in the XDMF order), the elset id
  !> and the phase of each element, gathered from every rank and ordered by global id.
  !> Collective.
  !>
  !> Arguments:
  !>   sim   (in) type(sim_writer_t)  the writer (HDF5 file open on the root)
  !>   mesh  (in) type(mesh_t)        rank-local mesh
  !> Locals:
  !>   n_own      integer(ik)         owned nodes
  !>   nen        integer(ik)         nodes per element
  !>   ids, eids  integer(ik)(:)      gathered global node and element ids (root)
  !>   vals, conn, elsets, phases  real(rk)(:)  gathered coordinates, element nodes, elset
  !>                                  ids, phases (root; integers as exact reals)
  !>   coo        real(rk)(:, :)      (3, nodes) coordinates by global id (root)
  !>   elts       real(rk)(:, :)      (nen, elts) element nodes by global element id (root)
  !>   el, ph     real(rk)(:, :)      (1, elts) elset ids, phases (root)
  !>   order      integer(ik)(:)      the XDMF node order
  subroutine h5_write_mesh(sim, mesh)

    type(sim_writer_t), intent(in) :: sim
    type(mesh_t), intent(in) :: mesh

    integer(ik) :: n_own, nen
    integer(ik), allocatable :: ids(:), eids(:), order(:)
    real(rk), allocatable :: vals(:), conn(:), elsets(:), phases(:)
    real(rk), allocatable :: coo(:, :), elts(:, :), el(:, :), ph(:, :)

    ! Every rank's owned nodes and elements on the root.
    n_own = mesh%num_nodes_owned
    nen = mesh%num_elt_nodes
    call par_gatherv(mesh%node_id_glob(1:n_own), ids)
    call par_gatherv(reshape(mesh%coo_s(:, 1:n_own), [3*n_own]), vals)
    call par_gatherv(mesh%elt_id_glob, eids)
    call par_gatherv(real(reshape(mesh%elt_nodes_glob, [nen*mesh%num_elts]), rk), conn)
    call par_gatherv(real(mesh%elset_id(mesh%elt_elset), rk), elsets)
    call par_gatherv(real(mesh%elt_phase, rk), phases)
    if (my_rank /= root_rank) return

    ! Placed by global id, then written.
    allocate (coo(3, mesh%num_nodes_glob), elts(nen, mesh%num_elts_glob))
    allocate (el(1, mesh%num_elts_glob), ph(1, mesh%num_elts_glob))
    call place(ids, vals, coo, 'sim_open', 'mesh coordinates')
    call place(eids, conn, elts, 'sim_open', 'mesh elements')
    call place(eids, elsets, el, 'sim_open', 'mesh elsets')
    call place(eids, phases, ph, 'sim_open', 'mesh phases')
    order = sim_xdmf_node_order(mesh%elt_type)
    call h5_write_real(sim%h5, '/mesh/coo', coo)
    call h5_write_int(sim%h5, '/mesh/elts', nint(elts, ik))
    call h5_write_int(sim%h5, '/mesh/elts_xdmf', nint(elts(order, :), ik) - 1_ik)
    call h5_write_int(sim%h5, '/mesh/elset', nint(el, ik))
    call h5_write_int(sim%h5, '/mesh/phase', nint(ph, ik))
    call h5_flush(sim%h5)

  end subroutine h5_write_mesh

  !> The XDMF (VTK) order of an element's nodes (step 9.1): node i of the XDMF element is node
  !> order(i) of FEPX's order (tet10_xdmf_order, hex20_xdmf_order).
  !>
  !> Arguments:
  !>   elt_type  (in) integer(ik)    Gmsh type of the element (11 or 17)
  !>   order     (out, result) integer(ik)(:)  the order (10 or 20 entries)
  function sim_xdmf_node_order(elt_type) result(order)

    integer(ik), intent(in) :: elt_type
    integer(ik), allocatable :: order(:)

    select case (elt_type)
    case (11_ik)
      order = tet10_xdmf_order
    case (17_ik)
      order = hex20_xdmf_order
    case default
      call internal_error('sim_xdmf_node_order', 'no XDMF node order for element type ' // &
        & int_text(int(elt_type)))
    end select

  end function sim_xdmf_node_order

  !> Write the XDMF file (step 9.1; XDMF 3, the XDMF Model and Format page, xdmf.org): a
  !> temporal collection with one grid per step written, each referring to the HDF5 datasets
  !> of its step (the datasets present: a variable is listed at the steps where it was
  !> written). Replaces the previous file. Root only.
  !>
  !> Arguments:
  !>   sim  (in) type(sim_writer_t)  the writer (HDF5 file open)
  !> Locals:
  !>   unit       integer            file unit
  !>   j, i       integer            step, variable
  !>   st         character(:), allocatable  the step number as text
  !>   path, var  character(:), allocatable  a dataset, a variable
  !>   rows, cols integer            its dimensions
  !>   kind       character(:), allocatable  its XDMF attribute type
  !>   topology   character(6)       Tet_10 or Hex_20
  subroutine write_xdmf(sim)

    type(sim_writer_t), intent(in) :: sim

    integer :: unit, j, i, rows, cols
    character(len=:), allocatable :: st, path, var, kind
    character(len=6) :: topology

    topology = 'Tet_10'
    if (sim%elt_type == 17_ik) topology = 'Hex_20'
    call delete_if_exists(sim%dir // '/' // xdmf_base)
    unit = open_new(sim%dir // '/' // xdmf_base)
    write (unit, '(a)') '<?xml version="1.0" ?>'
    write (unit, '(a)') '<Xdmf Version="3.0">'
    write (unit, '(a)') '  <Domain>'
    write (unit, '(a)') '    <Grid Name="simulation" GridType="Collection" ' // &
      & 'CollectionType="Temporal">'

    ! One grid per step written: time, topology, geometry, then the attributes.
    do j = 1, size(sim%h5_steps)
      st = int_text(sim%h5_steps(j))
      write (unit, '(a)') '      <Grid Name="step' // st // '" GridType="Uniform">'
      write (unit, '(a)') '        <Time Value="' // real_text(sim%h5_times(j)) // '"/>'
      write (unit, '(a)') '        <Topology TopologyType="' // trim(topology) // &
        & '" NumberOfElements="' // int_text(int(sim%general(2))) // '">'
      call data_item('/mesh/elts_xdmf', 'Int', '4')
      write (unit, '(a)') '        </Topology>'
      write (unit, '(a)') '        <Geometry GeometryType="XYZ">'
      path = '/results/nodes/coo/coo.step' // st
      if (.not. h5_exists(sim%h5, path)) path = '/mesh/coo'
      call data_item(path, 'Float', '8')
      write (unit, '(a)') '        </Geometry>'
      call attribute('elset', 'Scalar', 'Cell', '/mesh/elset', 'Int', '4')
      call attribute('phase', 'Scalar', 'Cell', '/mesh/phase', 'Int', '4')

      ! Node variables (vectors), then element variables.
      do i = 1, size(sim%node_vars)
        var = trim(sim%node_vars(i))
        path = '/results/nodes/' // var // '/' // var // '.step' // st
        if (h5_exists(sim%h5, path)) call attribute(var, 'Vector', 'Node', path, 'Float', '8')
      end do
      do i = 1, size(sim%elt_vars)
        var = trim(sim%elt_vars(i))
        if (var == 'forces' .or. var == 'convergence') cycle
        path = '/results/elts/' // var // '/' // var // '.step' // st
        if (.not. h5_exists(sim%h5, path)) cycle
        call h5_dims(sim%h5, path, rows, cols)
        kind = 'Matrix'
        if (cols == 1) kind = 'Scalar'
        if (cols == 3 .and. index(xdmf_vector_names, ' ' // var // ' ') > 0) kind = 'Vector'
        call attribute(var, kind, 'Cell', path, 'Float', '8')
      end do
      write (unit, '(a)') '      </Grid>'
    end do
    write (unit, '(a)') '    </Grid>'
    write (unit, '(a)') '  </Domain>'
    write (unit, '(a)') '</Xdmf>'
    close (unit)

  contains

    !> One attribute: its name, type and centre, and its data item.
    !>
    !> Arguments:
    !>   name, attr_type, center  (in) character  the XDMF attribute
    !>   item_path, number_type, precision  (in) character  its data item (data_item)
    subroutine attribute(name, attr_type, center, item_path, number_type, precision)

      character(len=*), intent(in) :: name, attr_type, center, item_path, number_type, &
        & precision

      write (unit, '(a)') '        <Attribute Name="' // name // '" AttributeType="' // &
        & attr_type // '" Center="' // center // '">'
      call data_item(item_path, number_type, precision)
      write (unit, '(a)') '        </Attribute>'

    end subroutine attribute

    !> One HDF5 data item: its dimensions (rows, columns), number type and precision, and the
    !> dataset (file name relative to the XDMF file).
    !>
    !> Arguments:
    !>   item_path    (in) character  the dataset
    !>   number_type  (in) character  Float or Int
    !>   precision    (in) character  bytes: 8 or 4
    !> Locals:
    !>   n, k         integer         its dimensions
    subroutine data_item(item_path, number_type, precision)

      character(len=*), intent(in) :: item_path, number_type, precision

      integer :: n, k

      call h5_dims(sim%h5, item_path, n, k)
      write (unit, '(a)') '          <DataItem Dimensions="' // int_text(n) // ' ' // &
        & int_text(k) // '" NumberType="' // number_type // '" Precision="' // precision // &
        & '" Format="HDF">' // h5_base // ':' // item_path // '</DataItem>'

    end subroutine data_item

  end subroutine write_xdmf

  !> The names of a list, separated by blanks (an HDF5 attribute).
  !>
  !> Arguments:
  !>   names  (in) character(*)(:)  the names
  !>   text   (out, result) character(:), allocatable
  !> Locals:
  !>   i      integer               index
  function list_text(names) result(text)

    character(len=*), intent(in) :: names(:)
    character(len=:), allocatable :: text

    integer :: i

    text = ''
    do i = 1, size(names)
      if (i > 1) text = text // ' '
      text = text // trim(names(i))
    end do

  end function list_text

  !> Place gathered values by global id: values(k*(j-1)+1 : k*j) belong to ids(j); every id
  !> 1..size(table, 2) must appear exactly once.
  !>
  !> Arguments:
  !>   ids      (in)  integer(ik)(:)   global ids
  !>   values   (in)  real(rk)(:)      k values per id
  !>   table    (out) real(rk)(:, :)   (k, number of ids) values by id
  !>   routine, name  (in)  character  for messages
  !> Locals:
  !>   seen     logical(:)             id already placed
  !>   k, j     integer(ik)            values per id, index
  subroutine place(ids, values, table, routine, name)

    integer(ik), intent(in) :: ids(:)
    real(rk), intent(in) :: values(:)
    real(rk), intent(out) :: table(:, :)
    character(len=*), intent(in) :: routine, name

    logical, allocatable :: seen(:)
    integer(ik) :: k, j

    k = size(table, 1)
    if (size(ids) /= size(table, 2)) call internal_error(routine, name // &
      & ': gathered entity count differs from the mesh')
    allocate (seen(size(table, 2)))
    seen = .false.
    do j = 1_ik, size(ids, kind=ik)
      if (ids(j) < 1_ik .or. ids(j) > size(table, 2)) call internal_error(routine, name // &
        & ': global id out of range')
      if (seen(ids(j))) call internal_error(routine, name // ': global id gathered twice')
      seen(ids(j)) = .true.
      table(:, ids(j)) = values(k*(j - 1_ik) + 1_ik:k*j)
    end do

  end subroutine place

  !> Create a directory and its parents (mkdir -p; Fortran has no directory statement).
  !>
  !> Arguments:
  !>   path  (in)  character  directory
  !> Locals:
  !>   stat, cmdstat  integer  exit status of the command, status of its execution
  subroutine make_dir(path)

    character(len=*), intent(in) :: path

    integer :: stat, cmdstat

    stat = 0
    cmdstat = 0
    call execute_command_line('mkdir -p "' // path // '"', exitstat=stat, cmdstat=cmdstat)
    if (stat /= 0 .or. cmdstat /= 0) call fatal_error('sim_open', 'cannot create directory ' // &
      & path)

  end subroutine make_dir

  !> Copy a file byte for byte (stream access).
  !>
  !> Arguments:
  !>   source, target  (in)  character  file names
  !>   replace         (in, optional) logical  replace an existing target (default: an error)
  !> Locals:
  !>   in, out, ios    integer          units, status
  !>   size_bytes      integer(ik)      file size
  !>   bytes           character(:), allocatable  the contents
  subroutine copy_file(source, target, replace)

    character(len=*), intent(in) :: source, target
    logical, intent(in), optional :: replace

    integer :: in, out, ios
    integer(ik) :: size_bytes
    character(len=:), allocatable :: bytes

    open (newunit=in, file=source, access='stream', form='unformatted', status='old', &
      & action='read', iostat=ios)
    if (ios /= 0) call fatal_error('sim_open', 'cannot read input file ' // source)
    inquire (unit=in, size=size_bytes)
    allocate (character(len=size_bytes) :: bytes)
    if (size_bytes > 0_ik) read (in) bytes
    close (in)
    if (present(replace)) then
      if (replace) call delete_if_exists(target)
    end if
    open (newunit=out, file=target, access='stream', form='unformatted', status='new', &
      & action='write', iostat=ios)
    if (ios /= 0) call fatal_error('sim_open', 'cannot write ' // target)
    write (out) bytes
    close (out)

  end subroutine copy_file

  !> Open a new formatted file for writing (an existing file is an error).
  !>
  !> Arguments:
  !>   name  (in)  character  file name
  !> Locals:
  !>   ios   integer          status
  function open_new(name) result(unit)

    character(len=*), intent(in) :: name
    integer :: unit

    integer :: ios

    open (newunit=unit, file=name, status='new', action='write', form='formatted', &
      & iostat=ios)
    if (ios /= 0) call fatal_error('sim_writer', 'cannot create ' // name)

  end function open_new

  !> Open an existing formatted file, keep its first num_lines lines, remove the rest
  !> (ENDFILE after them) and reopen it for appending (step 8.4: a file continued from a
  !> restart point; fewer lines than num_lines is an error).
  !>
  !> Arguments:
  !>   name       (in)  character  file name
  !>   num_lines  (in)  integer    lines to keep
  !>   unit       (out, result) integer  the unit, positioned at the end
  !> Locals:
  !>   ios, i     integer          status, line
  function open_truncated(name, num_lines) result(unit)

    character(len=*), intent(in) :: name
    integer, intent(in) :: num_lines
    integer :: unit

    integer :: ios, i

    ! Skip the lines kept, then end the file there.
    open (newunit=unit, file=name, status='old', action='readwrite', form='formatted', &
      & access='sequential', iostat=ios)
    if (ios /= 0) call fatal_error('sim_reopen', 'restart on: cannot open ' // name)
    do i = 1, num_lines
      read (unit, '(a)', iostat=ios)
      if (ios /= 0) call fatal_error('sim_reopen', 'restart on: ' // name // ' is shorter ' // &
        & 'than at the restart point')
    end do
    endfile (unit)
    close (unit)

    ! Reopened at its end for the rows of this run.
    open (newunit=unit, file=name, status='old', action='write', form='formatted', &
      & position='append', iostat=ios)
    if (ios /= 0) call fatal_error('sim_reopen', 'restart on: cannot reopen ' // name)

  end function open_truncated

  !> Delete a file if it exists (open, then close with status 'delete').
  !>
  !> Arguments:
  !>   name  (in)  character  file name
  !> Locals:
  !>   exists  logical  whether it exists
  !>   unit, ios  integer  unit, status
  subroutine delete_if_exists(name)

    character(len=*), intent(in) :: name

    logical :: exists
    integer :: unit, ios

    inquire (file=name, exist=exists)
    if (.not. exists) return
    open (newunit=unit, file=name, status='old', iostat=ios)
    if (ios == 0) close (unit, status='delete', iostat=ios)
    if (ios /= 0) call fatal_error('sim_writer', 'cannot remove ' // name)

  end subroutine delete_if_exists

  !> The file name without its directories.
  !>
  !> Arguments:
  !>   path  (in)  character  a path
  pure function base_name(path) result(name)

    character(len=*), intent(in) :: path
    character(len=:), allocatable :: name

    name = path(index(path, '/', back=.true.) + 1:)

  end function base_name

  !> Text of a real in sim_real_format, without blanks.
  !>
  !> Arguments:
  !>   value  (in)  real(rk)  the real
  !> Locals:
  !>   buf    character(32)   formatting buffer
  function real_text(value) result(text)

    real(rk), intent(in) :: value
    character(len=:), allocatable :: text

    character(len=32) :: buf

    write (buf, '(' // sim_real_format // ')') value
    text = trim(adjustl(buf))

  end function real_text

  !> Decimal text of an integer.
  !>
  !> Arguments:
  !>   value  (in)  integer  the integer
  !> Locals:
  !>   buf    character(24)  formatting buffer
  pure function int_text(value) result(text)

    integer, intent(in) :: value
    character(len=:), allocatable :: text

    character(len=24) :: buf

    write (buf, '(i0)') value
    text = trim(buf)

  end function int_text

end module sim_writer_mod
