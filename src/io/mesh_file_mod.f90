! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Neper mesh file (.msh) parser: reads a whole mesh file, and optional auxiliary files, in one
!> scan into a whole-mesh structure (run on rank 0 only; the distribution to the ranks is done by
!> mesh_distribute_mod).
!>
!> Source : port of v2 read_input_msh.f90 / read_input_msh2.f90 (sections, their layout, the
!>          element and faset tags, the auxiliary-file sections of simulation.opt), restructured
!>          as one scan (ROADMAP D32) with strict checks; file format: Neper documentation
!>          doc/fileformat.rst, "Mesh File (.msh)" (Neper 5.0). Contains no numerical method.
!>          v2-fixed CHANGES.md C20-C22 were found while porting it.
!>
!> Policy:
!>   - Only Neper msh file version 2.4 is read ($MeshVersion 2.4; ROADMAP D107): a missing or other
!>     version is a fatal error. $MeshFormat must be "2.2 0 8" (ASCII, 8-byte reals).
!>   - Every known section is parsed strictly (exact token counts, ids in range, footer line);
!>     an unknown section is skipped with a warning (D34). Blank lines between sections are
!>     ignored. A known section given twice in the mesh file is an error.
!>   - Sections are stored as read; references between sections (elset ids, element ids, node
!>     ids) are resolved and checked by mesh_file_finalize, so the order of the sections does not
!>     matter.
!>   - 3D elements: second-order tetrahedra (Gmsh type 11, 10 nodes) or hexahedra (type 17, 20
!>     nodes), not mixed; 3 tags (elset, elset, partition; Neper); ids consecutive. Lower-
!>     dimensional elements (Gmsh types 15, 1, 8, 2, 3, 9, 16, 10) are skipped. Element nodes
!>     are kept in the Gmsh order of the file, which is FEPX's internal element node order
!>     (ROADMAP D35; tet10_mod). v2 renumbered them when reading; FEPX does not.
!>   - Elsets are numbered 1..num_elsets in increasing order of their ids (elset_id).
!>   - $Fasets face lines: Neper writes the nodes of each boundary face (a 6-node triangle or an
!>     8-node quadrangle of its 2D mesh) in Gmsh order (corners first, then edge nodes; msh file
!>     version 2.4, ROADMAP D107), and the reader keeps them. The faces are assumed to be
!>     oriented out of the body, as Neper writes them; neither checked nor corrected (ROADMAP
!>     D37; Neper's x0, y0, z0 faces of a periodic single-cell mesh point into the body, a Neper
!>     error reported upstream, dev/neper_findings.md, A9).
!>   - Ids in element-level sections ($ElementOrientations, $ElementCrss, ...) are the 3D element
!>     numbers 1..num_elts (as Neper writes them); ids in $Fasets are file element ids (Neper
!>     doc: "3D elements adjacent to the boundary"), mapped to 3D element numbers. Each face's
!>     nodes must be nodes of its 3D element (which also detects files with the old faset ids).
!>   - Auxiliary files (v2's simulation.opt, orientation and phase files) may contain the
!>     orientation, group and value sections only; a section read from an auxiliary file
!>     replaces the same section read before (as v2's opt, ori and phase files), with a log line.
!>     They may also contain $ElsetCrySym (Neper 5 writes it first in its orientation files,
!>     step 8.5), which must then be the mesh file's (an error otherwise; kept if the mesh file
!>     has none).
!>   - Element orientations, if given, are used instead of elset orientations (as v2).
!>   - Orientations are converted to rot_c2s (orientation_mod) in mesh_file_finalize; the
!>     descriptor and convention are kept for output.
!>   - Global sizes are checked against huge(ik) (ROADMAP D28): 3 num_nodes (global dofs) and
!>     num_elts times the nodes per element.
!>
!> Public entities:
!>   mesh_file_t          whole-mesh content of the files
!>   mesh_nset_t          one node set
!>   mesh_faset_t         one face set (surface elements, each with its adjacent 3D element)
!>   mesh_values_t        one per-elset or per-element value section ($ElsetCrss, ...)
!>   mesh_file_read       parse a mesh file (.msh)
!>   mesh_file_read_aux   parse an auxiliary file
!>   mesh_file_finalize   resolve and check references, convert orientations (after the reads)
!>   mesh_file_load       mesh_file_read + mesh_file_read_aux for each auxiliary file + finalize
!>   num_value_sections, value_section_names   the value sections ($ElsetCrss, ..., per entity)
!>   elt_type_tet10, elt_type_hex20            Gmsh codes of the 3D element types
!>   supported_mesh_version                    "2.4"
module mesh_file_mod

  use kinds_mod, only: rk, ik, i8
  use error_mod, only: fatal_error, internal_error
  use log_mod, only: log_info, log_warning, log_file
  use text_file_mod, only: text_file_t, text_open, text_close, text_next_line, &
    & text_require_line, text_error, count_tokens, get_token, parse_ints, parse_reals, &
    & parse_int_token
  use orientation_mod, only: ori_descriptor_from_name, ori_convention_from_name, &
    & ori_num_params, rot_c2s_from_descriptor, ori_max_params
  use tet10_mod, only: tet10_gmsh_type
  use hex20_mod, only: hex20_gmsh_type

  implicit none

  private

  public :: mesh_file_t, mesh_nset_t, mesh_faset_t, mesh_values_t
  public :: mesh_file_read, mesh_file_read_aux, mesh_file_finalize, mesh_file_load
  public :: num_value_sections, value_section_names, value_section_is_elset
  public :: elt_type_tet10, elt_type_hex20, supported_mesh_version, elt_type_name

  !> Gmsh element type codes of the supported 3D elements (defined by the element modules).
  integer(ik), parameter :: elt_type_tet10 = tet10_gmsh_type
  integer(ik), parameter :: elt_type_hex20 = hex20_gmsh_type

  !> The only mesh file version read (ROADMAP D33).
  character(len=*), parameter :: supported_mesh_version = '2.4'

  !> The value sections of simulation.opt: v2's (per elset, then per element), then the
  !> initial network dislocation density of Kocks-Mecking hardening (step 7.2, D60; per elset,
  !> then per element), then the loop number density and loop diameter of irradiation
  !> softening (step 7.3, D61; each per elset, then per element).
  integer, parameter :: num_value_sections = 12
  character(len=24), parameter :: value_section_names(num_value_sections) = [character(len=24) :: &
    & 'ElsetCrss', 'ElsetCrssSat', 'ElsetHardCoeff', 'ElementCrss', 'ElementCrssSat', &
    & 'ElementHardCoeff', 'ElsetRho0', 'ElementRho0', 'ElsetIrrLoopDensity', &
    & 'ElementIrrLoopDensity', 'ElsetIrrLoopDiameter', 'ElementIrrLoopDiameter']
  logical, parameter :: value_section_is_elset(num_value_sections) = &
    & [.true., .true., .true., .false., .false., .false., .true., .false., .true., .false., &
    & .true., .false.]

  !> One node set.
  !>
  !> Components:
  !>   label  character(:), allocatable  nset label (e.g. "x0")
  !>   nodes  integer(ik)(:)             node ids (1..num_nodes), in file order
  type :: mesh_nset_t
    character(len=:), allocatable :: label
    integer(ik), allocatable :: nodes(:)
  end type mesh_nset_t

  !> One face set.
  !>
  !> Components:
  !>   label           character(:), allocatable  faset label (e.g. "z1")
  !>   num_face_nodes  integer(ik)                nodes per face (6: tri6, 8: quad8)
  !>   elts            integer(ik)(:)             3D element number (1..num_elts) of each face
  !>   nodes           integer(ik)(:,:)           (num_face_nodes, num faces) node ids, Gmsh
  !>                                              order of the face element, as Neper
  !>                                              orients it (module policy)
  type :: mesh_faset_t
    character(len=:), allocatable :: label
    integer(ik) :: num_face_nodes = 0_ik
    integer(ik), allocatable :: elts(:)
    integer(ik), allocatable :: nodes(:, :)
  end type mesh_faset_t

  !> One per-elset or per-element value section, after mesh_file_finalize.
  !>
  !> Components:
  !>   defined        logical         whether the section was read
  !>   source         character(:)    file it was read from
  !>   num_values_hdr integer(ik)     second header value (in v2: values per entity after the
  !>                                  expansion by crystal type; checked by the material input)
  !>   num_values     integer(ik)(:)  number of values given for each entity
  !>   values         real(rk)(:,:)   (max given, num entities) values, by elset index or element
  !>                                  number [MPa for strengths]
  type :: mesh_values_t
    logical :: defined = .false.
    character(len=:), allocatable :: source
    integer(ik) :: num_values_hdr = 0_ik
    integer(ik), allocatable :: num_values(:)
    real(rk), allocatable :: values(:, :)
  end type mesh_values_t

  !> A section of "id value..." lines as read (resolved by mesh_file_finalize).
  !>
  !> Components:
  !>   defined  logical         whether it was read
  !>   source   character(:)    file it was read from
  !>   header   character(:)    header token after the count (descriptor, or value count)
  !>   ids      integer(ik)(:)  id at the start of each line
  !>   num      integer(ik)(:)  number of values on each line
  !>   vals     real(rk)(:,:)   (max values, lines) values
  type :: id_table_t
    logical :: defined = .false.
    character(len=:), allocatable :: source
    character(len=:), allocatable :: header
    integer(ik), allocatable :: ids(:)
    integer(ik), allocatable :: num(:)
    real(rk), allocatable :: vals(:, :)
  end type id_table_t

  !> Whole-mesh content.
  !>
  !> Components (sizes; units):
  !>   mesh_name       character(:)   mesh file name
  !>   version         character(:)   $MeshVersion ("2.4")
  !>   domain          character(:)   $Domain (e.g. "cube"), or empty
  !>   elset_crysym    character(:)   $ElsetCrySym (e.g. "triclinic"), or empty
  !>   num_nodes       integer(ik)    number of nodes
  !>   coo             real(rk)(3, num_nodes)  node coordinates, sample frame [mm]
  !>   elt_type        integer(ik)    Gmsh type of the 3D elements (11 or 17)
  !>   num_elt_nodes   integer(ik)    nodes per 3D element (10 or 20)
  !>   num_elts        integer(ik)    number of 3D elements
  !>   elt_id_first    integer(ik)    file id of the first 3D element
  !>   elt_nodes       integer(ik)(num_elt_nodes, num_elts)  node ids, Gmsh order (ROADMAP D35)
  !>   elt_elset       integer(ik)(num_elts)  elset index (1..num_elsets) of each element
  !>   elt_part        integer(ik)(num_elts)  partition (1..num_parts), or 0 if unpartitioned
  !>   num_elsets      integer(ik)    number of elsets
  !>   elset_id        integer(ik)(num_elsets)  elset ids (file), increasing
  !>   num_parts       integer(ik)    number of partitions (0: unpartitioned)
  !>   node_part       integer(ik)(num_nodes)  partition of each node, or 0
  !>   nsets           mesh_nset_t(:)   node sets
  !>   fasets          mesh_faset_t(:)  face sets
  !>   num_periodicity integer(ik)    number of periodicity relations
  !>   per_secondary, per_primary  integer(ik)(num_periodicity)  node ids
  !>   per_vect        integer(ik)(3, num_periodicity)  scaled vector primary -> secondary
  !>   has_ori         logical        whether orientations were given
  !>   ori_from_elts   logical        element orientations (.true.) or elset orientations
  !>   ori_descriptor, ori_convention  integer  descriptor and convention codes (orientation_mod)
  !>                                  of the orientations used (for output)
  !>   elt_rot_c2s     real(rk)(3, 3, num_elts)  crystal-to-sample rotation of each element
  !>   num_phases      integer(ik)    number of phases (1 if no $Groups)
  !>   elset_phase     integer(ik)(num_elsets)  phase of each elset
  !>   elt_phase       integer(ik)(num_elts)  phase of each element
  !>   values          mesh_values_t(num_value_sections)  value sections ($ElsetCrss, ...)
  !>   (private, as read) raw_* tables, seen_* flags, faset file ids, partition counts
  type :: mesh_file_t
    character(len=:), allocatable :: mesh_name
    character(len=:), allocatable :: version
    character(len=:), allocatable :: domain
    character(len=:), allocatable :: elset_crysym
    integer(ik) :: num_nodes = 0_ik
    real(rk), allocatable :: coo(:, :)
    integer(ik) :: elt_type = 0_ik
    integer(ik) :: num_elt_nodes = 0_ik
    integer(ik) :: num_elts = 0_ik
    integer(ik) :: elt_id_first = 0_ik
    integer(ik), allocatable :: elt_nodes(:, :)
    integer(ik), allocatable :: elt_elset(:)
    integer(ik), allocatable :: elt_part(:)
    integer(ik) :: num_elsets = 0_ik
    integer(ik), allocatable :: elset_id(:)
    integer(ik) :: num_parts = 0_ik
    integer(ik), allocatable :: node_part(:)
    type(mesh_nset_t), allocatable :: nsets(:)
    type(mesh_faset_t), allocatable :: fasets(:)
    integer(ik) :: num_periodicity = 0_ik
    integer(ik), allocatable :: per_secondary(:)
    integer(ik), allocatable :: per_primary(:)
    integer(ik), allocatable :: per_vect(:, :)
    logical :: has_ori = .false.
    logical :: ori_from_elts = .false.
    integer :: ori_descriptor = 0
    integer :: ori_convention = 0
    real(rk), allocatable :: elt_rot_c2s(:, :, :)
    integer(ik) :: num_phases = 1_ik
    integer(ik), allocatable :: elset_phase(:)
    integer(ik), allocatable :: elt_phase(:)
    type(mesh_values_t) :: values(num_value_sections)
    ! As read, resolved by mesh_file_finalize.
    integer(ik), allocatable, private :: elt_elset_id(:)
    integer(ik), allocatable, private :: raw_node_part(:, :)
    type(id_table_t), private :: raw_elset_ori
    type(id_table_t), private :: raw_elt_ori
    type(id_table_t), private :: raw_groups
    type(id_table_t), private :: raw_values(num_value_sections)
    logical, private :: seen_format = .false.
    logical, private :: seen_version = .false.
    logical, private :: seen_nodes = .false.
    logical, private :: seen_elements = .false.
    logical, private :: seen_nsets = .false.
    logical, private :: seen_fasets = .false.
    logical, private :: seen_node_parts = .false.
    logical, private :: seen_periodicity = .false.
    logical, private :: finalized = .false.
  end type mesh_file_t

contains

  ! ==========================================================================================
  ! Public routines
  ! ==========================================================================================

  !> Read a mesh file, its auxiliary files and finalize (see the module policy).
  !>
  !> Arguments:
  !>   mesh       (out)           type(mesh_file_t)  whole-mesh content
  !>   mesh_name  (in)            character          mesh file name
  !>   aux_names  (in, optional)  character(*)(:)    auxiliary file names, read in order
  !> Locals:
  !>   i_file  integer  auxiliary file index
  subroutine mesh_file_load(mesh, mesh_name, aux_names)

    type(mesh_file_t), intent(out) :: mesh
    character(len=*), intent(in) :: mesh_name
    character(len=*), intent(in), optional :: aux_names(:)

    integer :: i_file

    ! Mesh file, then each auxiliary file, then resolve.
    call mesh_file_read(mesh, mesh_name)
    if (present(aux_names)) then
      do i_file = 1, size(aux_names)
        call mesh_file_read_aux(mesh, trim(aux_names(i_file)))
      end do
    end if
    call mesh_file_finalize(mesh)

  end subroutine mesh_file_load

  !> Parse a mesh file (one scan). Call mesh_file_finalize after this and any auxiliary files.
  !>
  !> Arguments:
  !>   mesh       (inout)  type(mesh_file_t)  whole-mesh content (should be default-initialized)
  !>   mesh_name  (in)     character          mesh file name
  !> Locals:
  !>   file  type(text_file_t)  the open file
  subroutine mesh_file_read(mesh, mesh_name)

    type(mesh_file_t), intent(inout) :: mesh
    character(len=*), intent(in) :: mesh_name

    type(text_file_t) :: file

    ! Parse every section, between the file lines of the log (step 9.7).
    mesh%mesh_name = mesh_name
    call log_file('i', 'Parsing file `' // mesh_name // "'...")
    call text_open(file, mesh_name)
    call read_sections(mesh, file, .false.)
    call text_close(file)
    call log_file('i', 'Parsed file `' // mesh_name // "'.")

  end subroutine mesh_file_read

  !> Parse an auxiliary file (orientation, group and value sections; see the module policy).
  !>
  !> Arguments:
  !>   mesh      (inout)  type(mesh_file_t)  whole-mesh content, after mesh_file_read
  !>   aux_name  (in)     character          auxiliary file name
  !> Locals:
  !>   file  type(text_file_t)  the open file
  subroutine mesh_file_read_aux(mesh, aux_name)

    type(mesh_file_t), intent(inout) :: mesh
    character(len=*), intent(in) :: aux_name

    type(text_file_t) :: file

    ! Parse its sections, between the file lines of the log (step 9.7).
    call log_file('i', 'Parsing file `' // aux_name // "'...")
    call text_open(file, aux_name)
    call read_sections(mesh, file, .true.)
    call text_close(file)
    call log_file('i', 'Parsed file `' // aux_name // "'.")

  end subroutine mesh_file_read_aux

  !> Resolve and check the references between sections, convert the orientations and set the
  !> defaults (see the module policy). Call once, after all reads.
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_file_t)  whole-mesh content
  !> Locals:
  !>   name  character(:), allocatable  mesh file name, for messages
  subroutine mesh_file_finalize(mesh)

    type(mesh_file_t), intent(inout) :: mesh

    character(len=:), allocatable :: name

    name = 'mesh file'
    if (allocated(mesh%mesh_name)) name = mesh%mesh_name
    if (mesh%finalized) call internal_error('mesh_file_finalize', name // ': already finalized')

    ! Required sections.
    if (.not. mesh%seen_format) call fatal_error('mesh_file_finalize', name // &
      & ': no $MeshFormat section')
    if (.not. mesh%seen_version) call fatal_error('mesh_file_finalize', name // &
      & ': no $MeshVersion section: meshes written before Neper msh file version ' // &
      & supported_mesh_version // ' are not supported (re-save the mesh with a current Neper)')
    if (.not. mesh%seen_nodes) call fatal_error('mesh_file_finalize', name // ': no $Nodes section')
    if (.not. mesh%seen_elements) call fatal_error('mesh_file_finalize', name // &
      & ': no $Elements section')

    ! Elements, elsets, partitions, then everything that refers to them.
    call resolve_elements(mesh, name)
    call resolve_partitions(mesh, name)
    call resolve_sets(mesh, name)
    call resolve_orientations(mesh, name)
    call resolve_groups(mesh, name)
    call resolve_values(mesh, name)

    mesh%finalized = .true.

  end subroutine mesh_file_finalize

  ! ==========================================================================================
  ! Section dispatch
  ! ==========================================================================================

  !> Read all sections of an open file.
  !>
  !> Arguments:
  !>   mesh    (inout)  type(mesh_file_t)  whole-mesh content
  !>   file    (inout)  type(text_file_t)  the open file
  !>   is_aux  (in)     logical            auxiliary file (restricted set of sections)
  !> Locals:
  !>   line     character(:), allocatable  current line
  !>   eof      logical                    end of file reached
  !>   section  character(:), allocatable  section name (without the leading $)
  !>   i_value  integer                    value-section index
  subroutine read_sections(mesh, file, is_aux)

    type(mesh_file_t), intent(inout) :: mesh
    type(text_file_t), intent(inout) :: file
    logical, intent(in) :: is_aux

    character(len=:), allocatable :: line
    logical :: eof
    character(len=:), allocatable :: section
    integer :: i_value

    do
      ! Next non-blank line: a section header "$Name".
      call text_next_line(file, line, eof)
      if (eof) exit
      if (len_trim(line) == 0) cycle
      if (line(1:1) /= '$' .or. count_tokens(line) /= 1) then
        call text_error(file, 'expected a section header ($Name), found "' // line // '"')
      end if
      section = line(2:)
      if (len(section) >= 3) then
        if (section(1:3) == 'End') call text_error(file, 'unexpected ' // line)
      end if

      ! Sections allowed in auxiliary files.
      i_value = value_section_index(section)
      if (i_value > 0) then
        call read_id_table(mesh%raw_values(i_value), file, section, is_aux, .false.)
        cycle
      end if
      select case (section)
      case ('ElsetOrientations')
        call read_id_table(mesh%raw_elset_ori, file, section, is_aux, .true.)
        cycle
      case ('ElementOrientations')
        call read_id_table(mesh%raw_elt_ori, file, section, is_aux, .true.)
        cycle
      case ('Groups')
        call read_groups(mesh%raw_groups, file, is_aux)
        cycle
      case ('ElsetCrySym')
        ! Crystal symmetry label: stored from the mesh file; in an auxiliary file, the mesh
        ! file's value if it has one.
        call read_single_line(file, section, line)
        if (allocated(mesh%elset_crysym)) then
          if (is_aux .and. line /= mesh%elset_crysym) then
            call text_error(file, '$ElsetCrySym "' // line // '" differs from the mesh file''s "' &
              & // mesh%elset_crysym // '"')
          end if
        end if
        if (.not. allocated(mesh%elset_crysym)) mesh%elset_crysym = line
        cycle
      end select

      ! Mesh-file sections.
      if (is_aux) call text_error(file, 'section ' // line // ' is not allowed in an auxiliary ' &
        & // 'file (only orientations, groups and value sections)')
      select case (section)
      case ('MeshFormat')
        call read_format(mesh, file)
      case ('MeshVersion')
        call read_version(mesh, file)
      case ('Domain')
        call read_single_line(file, section, mesh%domain)
      case ('Topology')
        call read_single_line(file, section, line)
      case ('Comments')
        call skip_section(file, section)
      case ('PhysicalNames')
        call read_physical_names(file)
      case ('Nodes')
        call read_nodes(mesh, file)
      case ('Elements')
        call read_elements(mesh, file)
      case ('NSets')
        call read_nsets(mesh, file)
      case ('Fasets')
        call read_fasets(mesh, file)
      case ('NodePartitions')
        call read_node_partitions(mesh, file)
      case ('Periodicity')
        call read_periodicity(mesh, file)
      case default
        ! Unknown section: skip it (D34).
        call log_warning('Skipping unknown mesh section ' // trim(line) // ' of `' // &
          & file%name // "'.", 1)
        call skip_section(file, section)
      end select
    end do

  end subroutine read_sections

  ! ==========================================================================================
  ! Section readers (each reads from the line after the header to the footer included)
  ! ==========================================================================================

  !> $MeshFormat: "2.2 0 8" (Gmsh format 2.2, ASCII, 8-byte reals).
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_file_t)  whole-mesh content
  !>   file  (inout)  type(text_file_t)  the open file
  !> Locals:
  !>   line  character(:), allocatable  current line
  subroutine read_format(mesh, file)

    type(mesh_file_t), intent(inout) :: mesh
    type(text_file_t), intent(inout) :: file

    character(len=:), allocatable :: line

    call check_once(file, mesh%seen_format, 'MeshFormat')
    call text_require_line(file, line, 'the mesh format line')
    if (count_tokens(line) /= 3 .or. get_token(line, 1) /= '2.2' .or. &
      & get_token(line, 2) /= '0' .or. get_token(line, 3) /= '8') then
      call text_error(file, 'mesh format must be "2.2 0 8" (Gmsh 2.2, ASCII, 8-byte reals), ' &
        & // 'found "' // line // '"')
    end if
    call read_footer(file, 'MeshFormat')

  end subroutine read_format

  !> $MeshVersion: must be the supported version (D107).
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_file_t)  whole-mesh content
  !>   file  (inout)  type(text_file_t)  the open file
  !> Locals:
  !>   line  character(:), allocatable  current line
  subroutine read_version(mesh, file)

    type(mesh_file_t), intent(inout) :: mesh
    type(text_file_t), intent(inout) :: file

    character(len=:), allocatable :: line

    call check_once(file, mesh%seen_version, 'MeshVersion')
    call text_require_line(file, line, 'the mesh version')
    if (count_tokens(line) /= 1) call text_error(file, 'expected one mesh version, found "' &
      & // line // '"')
    mesh%version = trim(adjustl(line))
    if (mesh%version /= supported_mesh_version) then
      call text_error(file, 'mesh file version ' // mesh%version // ' is not supported: FEPX ' &
        & // 'reads Neper msh file version ' // supported_mesh_version // ' only (re-save the ' &
        & // 'mesh with a current Neper)')
    end if
    call read_footer(file, 'MeshVersion')

  end subroutine read_version

  !> $Nodes: count, then "id x y z" lines with ids 1..count in order.
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_file_t)  whole-mesh content
  !>   file  (inout)  type(text_file_t)  the open file
  !> Locals:
  !>   line    character(:), allocatable  current line
  !>   count   integer(ik)                number of nodes
  !>   i_node  integer(ik)                node index
  !>   values  real(rk)(4)                id and coordinates of a line
  subroutine read_nodes(mesh, file)

    type(mesh_file_t), intent(inout) :: mesh
    type(text_file_t), intent(inout) :: file

    character(len=:), allocatable :: line
    integer(ik) :: count
    integer(ik) :: i_node
    real(rk) :: values(4)

    call check_once(file, mesh%seen_nodes, 'Nodes')
    count = read_count(file, 'number of nodes')

    ! Global dofs must fit in ik (D28).
    call check_size(file, 3_i8*int(count, i8), 'number of degrees of freedom (3 x nodes)')
    mesh%num_nodes = count
    allocate (mesh%coo(3, count))

    ! One line per node; ids must be 1..count in order (Neper doc).
    do i_node = 1_ik, count
      call text_require_line(file, line, 'a node line')
      call parse_reals(file, line, values, 'node line')
      if (parse_int_token(file, get_token(line, 1), 'node id') /= i_node) then
        call text_error(file, 'node ids must be 1, 2, ... in order: expected ' // &
          & int_text(i_node) // ', found ' // get_token(line, 1))
      end if
      mesh%coo(:, i_node) = values(2:4)
    end do
    call read_footer(file, 'Nodes')

  end subroutine read_nodes

  !> $Elements: count, then "id type num_tags tags... nodes..." lines. 3D elements (11 or 17)
  !> are stored, lower-dimensional ones skipped (module policy).
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_file_t)  whole-mesh content
  !>   file  (inout)  type(text_file_t)  the open file
  !> Locals:
  !>   line        character(:), allocatable  current line
  !>   count       integer(ik)                number of element lines
  !>   i_line      integer(ik)                line index
  !>   num_tokens  integer                    tokens on a line
  !>   fields      integer(ik)(:)             the integers of a line
  !>   elt_type    integer(ik)                Gmsh type of the line's element
  !>   num_3d      integer(ik)                3D elements stored so far
  !>   nodes_3d    integer(ik)(:,:)           3D element nodes, Gmsh order
  !>   elset_3d    integer(ik)(:)             3D element elset ids
  !>   part_3d     integer(ik)(:)             3D element partitions
  subroutine read_elements(mesh, file)

    type(mesh_file_t), intent(inout) :: mesh
    type(text_file_t), intent(inout) :: file

    character(len=:), allocatable :: line
    integer(ik) :: count
    integer(ik) :: i_line
    integer :: num_tokens
    integer(ik), allocatable :: fields(:)
    integer(ik) :: elt_type
    integer(ik) :: num_3d
    integer(ik), allocatable :: nodes_3d(:, :)
    integer(ik), allocatable :: elset_3d(:)
    integer(ik), allocatable :: part_3d(:)

    call check_once(file, mesh%seen_elements, 'Elements')
    count = read_count(file, 'number of elements')
    num_3d = 0_ik
    allocate (nodes_3d(0, 0), elset_3d(0), part_3d(0))

    do i_line = 1_ik, count
      ! Read every integer of the line.
      call text_require_line(file, line, 'an element line')
      num_tokens = count_tokens(line)
      if (num_tokens < 3) call text_error(file, 'element line: expected id, type and number ' &
        & // 'of tags, found "' // line // '"')
      if (allocated(fields)) deallocate (fields)
      allocate (fields(num_tokens))
      call parse_ints(file, line, fields, 'element line')
      elt_type = fields(2)

      select case (elt_type)
      case (elt_type_tet10, elt_type_hex20)
        ! A 3D element: the type must not change, 3 tags, then its nodes.
        if (num_3d == 0_ik) then
          mesh%elt_type = elt_type
          mesh%num_elt_nodes = merge(10_ik, 20_ik, elt_type == elt_type_tet10)
          mesh%elt_id_first = fields(1)
          call check_size(file, int(count, i8)*int(mesh%num_elt_nodes, i8), &
            & 'element connectivity size (elements x nodes per element)')
          deallocate (nodes_3d, elset_3d, part_3d)
          allocate (nodes_3d(mesh%num_elt_nodes, count), elset_3d(count), part_3d(count))
        else if (elt_type /= mesh%elt_type) then
          call text_error(file, 'mixed 3D element types ' // int_text(mesh%elt_type) // &
            & ' and ' // int_text(elt_type) // ' are not supported')
        end if
        if (fields(3) /= 3_ik .or. num_tokens /= 6 + mesh%num_elt_nodes) then
          call text_error(file, 'element line: expected 3 tags (elset, elset, partition) ' &
            & // 'and ' // int_text(mesh%num_elt_nodes) // ' nodes, found "' // line // '"')
        end if
        num_3d = num_3d + 1_ik
        if (fields(1) /= mesh%elt_id_first + num_3d - 1_ik) then
          call text_error(file, '3D element ids must be consecutive: expected ' // &
            & int_text(mesh%elt_id_first + num_3d - 1_ik) // ', found ' // int_text(fields(1)))
        end if
        elset_3d(num_3d) = fields(4)
        part_3d(num_3d) = fields(6)
        ! Nodes in the file's (Gmsh) order, the internal element node order (D35).
        nodes_3d(:, num_3d) = fields(7:)
      case (15, 1, 8, 2, 3, 9, 16, 10)
        ! Lower-dimensional element: skipped.
        continue
      case default
        call text_error(file, 'element type ' // int_text(elt_type) // ' is not supported ' // &
          & '(3D elements must be second-order tetrahedra, type 11, or hexahedra, type 17)')
      end select
    end do
    call read_footer(file, 'Elements')

    ! At least one 3D element.
    if (num_3d == 0_ik) call text_error(file, 'no second-order 3D elements (type 11 or 17) ' &
      & // 'in $Elements')

    ! Keep the 3D elements only.
    mesh%num_elts = num_3d
    mesh%elt_nodes = nodes_3d(:, 1:num_3d)
    mesh%elt_elset_id = elset_3d(1:num_3d)
    mesh%elt_part = part_3d(1:num_3d)

  end subroutine read_elements

  !> $NSets: count, then for each nset its label, node count and node ids (one per line).
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_file_t)  whole-mesh content
  !>   file  (inout)  type(text_file_t)  the open file
  !> Locals:
  !>   line       character(:), allocatable  current line
  !>   num_nsets  integer(ik)                number of nsets
  !>   i_nset     integer(ik)                nset index
  !>   num        integer(ik)                nodes of an nset
  !>   i_node     integer(ik)                node index within the nset
  !>   id         integer(ik)(1)             a node id
  subroutine read_nsets(mesh, file)

    type(mesh_file_t), intent(inout) :: mesh
    type(text_file_t), intent(inout) :: file

    character(len=:), allocatable :: line
    integer(ik) :: num_nsets
    integer(ik) :: i_nset
    integer(ik) :: num
    integer(ik) :: i_node
    integer(ik) :: id(1)

    call check_once(file, mesh%seen_nsets, 'NSets')
    num_nsets = read_count(file, 'number of nsets')
    allocate (mesh%nsets(num_nsets))
    do i_nset = 1_ik, num_nsets
      mesh%nsets(i_nset)%label = read_label(file, 'nset label')
      num = read_count(file, 'number of nodes of nset ' // mesh%nsets(i_nset)%label)
      allocate (mesh%nsets(i_nset)%nodes(num))
      do i_node = 1_ik, num
        call text_require_line(file, line, 'an nset node id')
        call parse_ints(file, line, id, 'nset node id')
        mesh%nsets(i_nset)%nodes(i_node) = id(1)
      end do
    end do
    call read_footer(file, 'NSets')

  end subroutine read_nsets

  !> $Fasets: count, then for each faset its label, face count and "elt_id node..." lines.
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_file_t)  whole-mesh content
  !>   file  (inout)  type(text_file_t)  the open file
  !> Locals:
  !>   line        character(:), allocatable  current line
  !>   num_fasets  integer(ik)                number of fasets
  !>   i_faset     integer(ik)                faset index
  !>   num         integer(ik)                faces of a faset
  !>   i_face      integer(ik)                face index
  !>   num_nodes   integer(ik)                nodes per face (from the first face line)
  !>   fields      integer(ik)(:)             the integers of a face line
  subroutine read_fasets(mesh, file)

    type(mesh_file_t), intent(inout) :: mesh
    type(text_file_t), intent(inout) :: file

    character(len=:), allocatable :: line
    integer(ik) :: num_fasets
    integer(ik) :: i_faset
    integer(ik) :: num
    integer(ik) :: i_face
    integer(ik) :: num_nodes
    integer(ik), allocatable :: fields(:)

    call check_once(file, mesh%seen_fasets, 'Fasets')
    num_fasets = read_count(file, 'number of fasets')
    allocate (mesh%fasets(num_fasets))
    do i_faset = 1_ik, num_fasets
      mesh%fasets(i_faset)%label = read_label(file, 'faset label')
      num = read_count(file, 'number of faces of faset ' // mesh%fasets(i_faset)%label)
      num_nodes = 0_ik
      allocate (mesh%fasets(i_faset)%elts(num))
      do i_face = 1_ik, num
        call text_require_line(file, line, 'a faset face line')
        ! The first line sets the nodes per face (6: tri6, 8: quad8).
        if (i_face == 1_ik) then
          num_nodes = int(count_tokens(line), ik) - 1_ik
          if (num_nodes /= 6_ik .and. num_nodes /= 8_ik) then
            call text_error(file, 'faset face line: expected an element id and 6 or 8 nodes, ' &
              & // 'found "' // line // '"')
          end if
          mesh%fasets(i_faset)%num_face_nodes = num_nodes
          allocate (mesh%fasets(i_faset)%nodes(num_nodes, num))
          allocate (fields(num_nodes + 1_ik))
        end if
        call parse_ints(file, line, fields, 'faset face line')
        mesh%fasets(i_faset)%elts(i_face) = fields(1)
        ! The nodes of the face, in the Gmsh order of the face element (module policy).
        mesh%fasets(i_faset)%nodes(:, i_face) = fields(2_ik:num_nodes + 1_ik)
      end do
      ! An empty faset (Neper writes them for periodic meshes) has no face line: no nodes,
      ! a (0, 0) table, so that it is allocated like any other.
      if (.not. allocated(mesh%fasets(i_faset)%nodes)) &
        & allocate (mesh%fasets(i_faset)%nodes(0, 0))
      if (allocated(fields)) deallocate (fields)
    end do
    call read_footer(file, 'Fasets')

  end subroutine read_fasets

  !> $NodePartitions: count, then "node_id partition" lines.
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_file_t)  whole-mesh content
  !>   file  (inout)  type(text_file_t)  the open file
  !> Locals:
  !>   line    character(:), allocatable  current line
  !>   count   integer(ik)                number of lines
  !>   i_line  integer(ik)                line index
  !>   fields  integer(ik)(2)             node id and partition
  subroutine read_node_partitions(mesh, file)

    type(mesh_file_t), intent(inout) :: mesh
    type(text_file_t), intent(inout) :: file

    character(len=:), allocatable :: line
    integer(ik) :: count
    integer(ik) :: i_line
    integer(ik) :: fields(2)

    call check_once(file, mesh%seen_node_parts, 'NodePartitions')
    count = read_count(file, 'number of node partition lines')
    allocate (mesh%raw_node_part(2, count))
    do i_line = 1_ik, count
      call text_require_line(file, line, 'a node partition line')
      call parse_ints(file, line, fields, 'node partition line')
      mesh%raw_node_part(:, i_line) = fields
    end do
    call read_footer(file, 'NodePartitions')

  end subroutine read_node_partitions

  !> $Periodicity: count, then "secondary primary px py pz" lines.
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_file_t)  whole-mesh content
  !>   file  (inout)  type(text_file_t)  the open file
  !> Locals:
  !>   line    character(:), allocatable  current line
  !>   count   integer(ik)                number of relations
  !>   i_line  integer(ik)                line index
  !>   fields  integer(ik)(5)             the integers of a line
  subroutine read_periodicity(mesh, file)

    type(mesh_file_t), intent(inout) :: mesh
    type(text_file_t), intent(inout) :: file

    character(len=:), allocatable :: line
    integer(ik) :: count
    integer(ik) :: i_line
    integer(ik) :: fields(5)

    call check_once(file, mesh%seen_periodicity, 'Periodicity')
    count = read_count(file, 'number of periodicity relations')
    mesh%num_periodicity = count
    allocate (mesh%per_secondary(count), mesh%per_primary(count), mesh%per_vect(3, count))
    do i_line = 1_ik, count
      call text_require_line(file, line, 'a periodicity line')
      call parse_ints(file, line, fields, 'periodicity line')
      if (any(abs(fields(3:5)) > 1_ik) .or. all(fields(3:5) == 0_ik)) then
        call text_error(file, 'periodicity vector components must be -1, 0 or 1, not all 0')
      end if
      mesh%per_secondary(i_line) = fields(1)
      mesh%per_primary(i_line) = fields(2)
      mesh%per_vect(:, i_line) = fields(3:5)
    end do
    call read_footer(file, 'Periodicity')

  end subroutine read_periodicity

  !> $PhysicalNames: count, then one line per name (not used; Gmsh information).
  !>
  !> Arguments:
  !>   file  (inout)  type(text_file_t)  the open file
  !> Locals:
  !>   line    character(:), allocatable  current line
  !>   count   integer(ik)                number of names
  !>   i_line  integer(ik)                line index
  subroutine read_physical_names(file)

    type(text_file_t), intent(inout) :: file

    character(len=:), allocatable :: line
    integer(ik) :: count
    integer(ik) :: i_line

    count = read_count(file, 'number of physical names')
    do i_line = 1_ik, count
      call text_require_line(file, line, 'a physical name line')
      if (count_tokens(line) /= 3) call text_error(file, 'physical name line: expected ' // &
        & 'dimension, id and name, found "' // line // '"')
    end do
    call read_footer(file, 'PhysicalNames')

  end subroutine read_physical_names

  !> A section of "id value..." lines with a header "count token": orientations ("count
  !> descriptor:convention", lines "id params...") and value sections ("count num_values",
  !> lines "id values..."; the number of values may vary between lines, as in v2's shorthands).
  !>
  !> Arguments:
  !>   table     (inout)  type(id_table_t)   storage of the section
  !>   file      (inout)  type(text_file_t)  the open file
  !>   section   (in)     character          section name (without $)
  !>   is_aux    (in)     logical            read from an auxiliary file (replaces a definition)
  !>   is_ori    (in)     logical            orientation section (header token: descriptor)
  !> Locals:
  !>   line        character(:), allocatable  current line
  !>   count       integer(ik)                number of lines
  !>   i_line      integer(ik)                line index
  !>   num_tokens  integer                    tokens on a line
  !>   values      real(rk)(:)                the numbers of a line (id first)
  !>   max_values  integer(ik)                values per line, allocated width
  subroutine read_id_table(table, file, section, is_aux, is_ori)

    type(id_table_t), intent(inout) :: table
    type(text_file_t), intent(inout) :: file
    character(len=*), intent(in) :: section
    logical, intent(in) :: is_aux
    logical, intent(in) :: is_ori

    character(len=:), allocatable :: line
    integer(ik) :: count
    integer(ik) :: i_line
    integer :: num_tokens
    real(rk), allocatable :: values(:)
    integer(ik) :: max_values

    ! A second definition: an error in the mesh file, a replacement from an auxiliary file.
    call check_redefinition(table, file, section, is_aux)

    ! Header: count and token.
    call text_require_line(file, line, 'the $' // section // ' header')
    if (count_tokens(line) /= 2) then
      if (is_ori) then
        call text_error(file, '$' // section // ' header: expected count and ' // &
          & 'descriptor:convention, found "' // line // '"')
      else
        call text_error(file, '$' // section // ' header: expected count and number of ' // &
          & 'values, found "' // line // '"')
      end if
    end if
    count = parse_int_token(file, get_token(line, 1), '$' // section // ' count')
    if (count < 0_ik) call text_error(file, 'negative count')
    table%header = get_token(line, 2)
    if (.not. is_ori) then
      if (parse_int_token(file, table%header, 'number of values') < 1_ik) then
        call text_error(file, 'number of values must be positive')
      end if
    end if

    ! Lines: an id and at least one value; widths up to ori_max_params or 64 values.
    max_values = merge(int(ori_max_params, ik), 64_ik, is_ori)
    table%defined = .true.
    table%source = file%name
    if (allocated(table%ids)) deallocate (table%ids, table%num, table%vals)
    allocate (table%ids(count), table%num(count), table%vals(max_values, count))
    table%vals = 0.0_rk
    do i_line = 1_ik, count
      call text_require_line(file, line, 'a $' // section // ' line')
      num_tokens = count_tokens(line)
      if (num_tokens < 2 .or. num_tokens - 1 > max_values) then
        call text_error(file, '$' // section // ' line: expected an id and 1 to ' // &
          & int_text(max_values) // ' values, found "' // line // '"')
      end if
      if (allocated(values)) deallocate (values)
      allocate (values(num_tokens))
      call parse_reals(file, line, values, '$' // section // ' line')
      table%ids(i_line) = parse_int_token(file, get_token(line, 1), '$' // section // ' id')
      table%num(i_line) = int(num_tokens - 1, ik)
      table%vals(1:num_tokens - 1, i_line) = values(2:)
    end do
    call read_footer(file, section)

  end subroutine read_id_table

  !> $Groups: "elset", count, then "elset_id group" lines.
  !>
  !> Arguments:
  !>   table   (inout)  type(id_table_t)   storage (ids: elset ids; vals(1,:): groups)
  !>   file    (inout)  type(text_file_t)  the open file
  !>   is_aux  (in)     logical            read from an auxiliary file (replaces a definition)
  !> Locals:
  !>   line    character(:), allocatable  current line
  !>   count   integer(ik)                number of lines
  !>   i_line  integer(ik)                line index
  !>   fields  integer(ik)(2)             elset id and group
  subroutine read_groups(table, file, is_aux)

    type(id_table_t), intent(inout) :: table
    type(text_file_t), intent(inout) :: file
    logical, intent(in) :: is_aux

    character(len=:), allocatable :: line
    integer(ik) :: count
    integer(ik) :: i_line
    integer(ik) :: fields(2)

    call check_redefinition(table, file, 'Groups', is_aux)

    ! The grouped entity must be elsets (Neper doc).
    call text_require_line(file, line, 'the group entity')
    if (trim(adjustl(line)) /= 'elset') call text_error(file, '$Groups: the group entity must ' &
      & // 'be "elset", found "' // line // '"')
    count = read_count(file, 'number of groups lines')

    table%defined = .true.
    table%source = file%name
    table%header = 'elset'
    if (allocated(table%ids)) deallocate (table%ids, table%num, table%vals)
    allocate (table%ids(count), table%num(count), table%vals(1, count))
    do i_line = 1_ik, count
      call text_require_line(file, line, 'a groups line')
      call parse_ints(file, line, fields, 'groups line')
      table%ids(i_line) = fields(1)
      table%num(i_line) = 1_ik
      table%vals(1, i_line) = real(fields(2), rk)
    end do
    call read_footer(file, 'Groups')

  end subroutine read_groups

  !> A section of one line (e.g. $Domain), returned.
  !>
  !> Arguments:
  !>   file     (inout)  type(text_file_t)          the open file
  !>   section  (in)     character                  section name (without $)
  !>   value    (out)    character(:), allocatable  the line
  subroutine read_single_line(file, section, value)

    type(text_file_t), intent(inout) :: file
    character(len=*), intent(in) :: section
    character(len=:), allocatable, intent(out) :: value

    call text_require_line(file, value, 'the $' // section // ' value')
    value = trim(adjustl(value))
    call read_footer(file, section)

  end subroutine read_single_line

  !> Skip a section up to and including its footer "$End<section>".
  !>
  !> Arguments:
  !>   file     (inout)  type(text_file_t)  the open file
  !>   section  (in)     character          section name (without $)
  !> Locals:
  !>   line  character(:), allocatable  current line
  subroutine skip_section(file, section)

    type(text_file_t), intent(inout) :: file
    character(len=*), intent(in) :: section

    character(len=:), allocatable :: line

    do
      call text_require_line(file, line, '$End' // section)
      if (line == '$End' // section) exit
    end do

  end subroutine skip_section

  ! ==========================================================================================
  ! Resolution (mesh_file_finalize)
  ! ==========================================================================================

  !> Elements: node ids in range, every node used, elsets numbered by increasing id.
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_file_t)  whole-mesh content
  !>   name  (in)     character          mesh file name, for messages
  !> Locals:
  !>   used        logical(:)       whether each node belongs to a 3D element
  !>   i_elt       integer(ik)      element index
  !>   max_id      integer(ik)      largest elset id
  !>   id_to_idx   integer(ik)(:)   elset index of each elset id (0 if unused)
  !>   i_id        integer(ik)      elset id
  subroutine resolve_elements(mesh, name)

    type(mesh_file_t), intent(inout) :: mesh
    character(len=*), intent(in) :: name

    logical, allocatable :: used(:)
    integer(ik) :: i_elt
    integer(ik) :: max_id
    integer(ik), allocatable :: id_to_idx(:)
    integer(ik) :: i_id

    ! Node ids of the elements in range; every node belongs to a 3D element (an unused node
    ! would make the stiffness matrix singular).
    if (any(mesh%elt_nodes < 1_ik) .or. any(mesh%elt_nodes > mesh%num_nodes)) then
      call fatal_error('mesh_file_finalize', name // ': $Elements refers to a node id outside ' &
        & // '1..' // int_text(mesh%num_nodes))
    end if
    allocate (used(mesh%num_nodes))
    used = .false.
    do i_elt = 1_ik, mesh%num_elts
      used(mesh%elt_nodes(:, i_elt)) = .true.
    end do
    if (.not. all(used)) then
      call fatal_error('mesh_file_finalize', name // ': node ' // &
        & int_text(int(findloc(used, .false., dim=1), ik)) // ' belongs to no 3D element')
    end if

    ! Elsets: ids must be positive; numbered 1..num_elsets in increasing id order.
    if (any(mesh%elt_elset_id < 1_ik)) call fatal_error('mesh_file_finalize', name // &
      & ': 3D element with a non-positive elset id')
    max_id = maxval(mesh%elt_elset_id)
    allocate (id_to_idx(max_id))
    id_to_idx = 0_ik
    id_to_idx(mesh%elt_elset_id) = 1_ik
    mesh%num_elsets = int(count(id_to_idx > 0_ik), ik)
    allocate (mesh%elset_id(mesh%num_elsets))
    mesh%num_elsets = 0_ik
    do i_id = 1_ik, max_id
      if (id_to_idx(i_id) > 0_ik) then
        mesh%num_elsets = mesh%num_elsets + 1_ik
        mesh%elset_id(mesh%num_elsets) = i_id
        id_to_idx(i_id) = mesh%num_elsets
      end if
    end do
    mesh%elt_elset = id_to_idx(mesh%elt_elset_id)

  end subroutine resolve_elements

  !> Partitions: all 3D elements partitioned (1..num_parts) or none (0); node partitions given
  !> for every node exactly when the elements are partitioned.
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_file_t)  whole-mesh content
  !>   name  (in)     character          mesh file name, for messages
  !> Locals:
  !>   i_line  integer(ik)  node partition line
  !>   seen    logical(:)   whether each node has a partition line
  subroutine resolve_partitions(mesh, name)

    type(mesh_file_t), intent(inout) :: mesh
    character(len=*), intent(in) :: name

    integer(ik) :: i_line
    logical, allocatable :: seen(:)

    allocate (mesh%node_part(mesh%num_nodes))
    mesh%node_part = 0_ik

    ! Unpartitioned mesh: every partition tag 0, and no $NodePartitions.
    if (all(mesh%elt_part == 0_ik)) then
      mesh%num_parts = 0_ik
      if (mesh%seen_node_parts) call fatal_error('mesh_file_finalize', name // &
        & ': $NodePartitions given, but the 3D elements are not partitioned')
      return
    end if

    ! Partitioned mesh: element partitions 1..num_parts.
    if (any(mesh%elt_part < 1_ik)) call fatal_error('mesh_file_finalize', name // &
      & ': some 3D elements have a partition and others not')
    mesh%num_parts = maxval(mesh%elt_part)

    ! Node partitions: one line per node, partitions in 1..num_parts.
    if (.not. mesh%seen_node_parts) call fatal_error('mesh_file_finalize', name // &
      & ': the 3D elements are partitioned, but there is no $NodePartitions section')
    if (size(mesh%raw_node_part, 2) /= mesh%num_nodes) call fatal_error('mesh_file_finalize', &
      & name // ': $NodePartitions has ' // int_text(int(size(mesh%raw_node_part, 2), ik)) // &
      & ' lines for ' // int_text(mesh%num_nodes) // ' nodes')
    allocate (seen(mesh%num_nodes))
    seen = .false.
    do i_line = 1_ik, mesh%num_nodes
      associate (id => mesh%raw_node_part(1, i_line), part => mesh%raw_node_part(2, i_line))
        if (id < 1_ik .or. id > mesh%num_nodes) call fatal_error('mesh_file_finalize', name // &
          & ': $NodePartitions: node id ' // int_text(id) // ' out of range')
        if (seen(id)) call fatal_error('mesh_file_finalize', name // ': $NodePartitions: node ' &
          & // int_text(id) // ' given twice')
        if (part < 1_ik .or. part > mesh%num_parts) call fatal_error('mesh_file_finalize', &
          & name // ': $NodePartitions: partition ' // int_text(part) // ' of node ' // &
          & int_text(id) // ' outside 1..' // int_text(mesh%num_parts))
        seen(id) = .true.
        mesh%node_part(id) = part
      end associate
    end do
    deallocate (mesh%raw_node_part)

  end subroutine resolve_partitions

  !> Node sets, face sets and periodicity: ids in range; faset element ids mapped to 3D element
  !> numbers, each face's nodes must be nodes of its element.
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_file_t)  whole-mesh content
  !>   name  (in)     character          mesh file name, for messages
  !> Locals:
  !>   i_set   integer(ik)  set index
  !>   i_face  integer(ik)  face index
  !>   i_elt   integer(ik)  3D element number
  !>   i_node  integer(ik)  face node index
  subroutine resolve_sets(mesh, name)

    type(mesh_file_t), intent(inout) :: mesh
    character(len=*), intent(in) :: name

    integer(ik) :: i_set
    integer(ik) :: i_face
    integer(ik) :: i_elt
    integer(ik) :: i_node

    ! Missing optional sections become empty.
    if (.not. allocated(mesh%nsets)) allocate (mesh%nsets(0))
    if (.not. allocated(mesh%fasets)) allocate (mesh%fasets(0))
    if (.not. allocated(mesh%per_secondary)) then
      allocate (mesh%per_secondary(0), mesh%per_primary(0), mesh%per_vect(3, 0))
    end if

    ! Node sets.
    do i_set = 1_ik, int(size(mesh%nsets), ik)
      if (any(mesh%nsets(i_set)%nodes < 1_ik) .or. &
        & any(mesh%nsets(i_set)%nodes > mesh%num_nodes)) then
        call fatal_error('mesh_file_finalize', name // ': nset ' // mesh%nsets(i_set)%label // &
          & ' has a node id outside 1..' // int_text(mesh%num_nodes))
      end if
    end do

    ! Face sets: file element ids -> 3D element numbers; face nodes among the element's nodes.
    do i_set = 1_ik, int(size(mesh%fasets), ik)
      associate (faset => mesh%fasets(i_set))
        if (size(faset%elts) > 0) then
          if (faset%num_face_nodes /= merge(6_ik, 8_ik, mesh%elt_type == elt_type_tet10)) then
            call fatal_error('mesh_file_finalize', name // ': faset ' // faset%label // ' has ' &
              & // int_text(faset%num_face_nodes) // '-node faces, which do not match the ' // &
              & elt_type_name(mesh%elt_type) // ' elements')
          end if
        end if
        do i_face = 1_ik, int(size(faset%elts), ik)
          i_elt = faset%elts(i_face) - mesh%elt_id_first + 1_ik
          if (i_elt < 1_ik .or. i_elt > mesh%num_elts) then
            call fatal_error('mesh_file_finalize', name // ': faset ' // faset%label // &
              & ' refers to element id ' // int_text(faset%elts(i_face)) // ', which is not ' &
              & // 'a 3D element')
          end if
          faset%elts(i_face) = i_elt
          do i_node = 1_ik, faset%num_face_nodes
            if (.not. any(mesh%elt_nodes(:, i_elt) == faset%nodes(i_node, i_face))) then
              call fatal_error('mesh_file_finalize', name // ': faset ' // faset%label // &
                & ': node ' // int_text(faset%nodes(i_node, i_face)) // ' of face ' // &
                & int_text(i_face) // ' is not a node of its element (id ' // &
                & int_text(i_elt + mesh%elt_id_first - 1_ik) // ')')
            end if
          end do
        end do
      end associate
    end do

    ! Periodicity: node ids in range.
    if (any(mesh%per_secondary < 1_ik) .or. any(mesh%per_secondary > mesh%num_nodes) .or. &
      & any(mesh%per_primary < 1_ik) .or. any(mesh%per_primary > mesh%num_nodes)) then
      call fatal_error('mesh_file_finalize', name // ': $Periodicity refers to a node id ' // &
        & 'outside 1..' // int_text(mesh%num_nodes))
    end if

  end subroutine resolve_sets

  !> Orientations: element orientations if given, else elset orientations, converted to
  !> rot_c2s for every element; descriptor and convention kept.
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_file_t)  whole-mesh content
  !>   name  (in)     character          mesh file name, for messages
  !> Locals:
  !>   order        integer(ik)(:)       entity index of each line
  !>   rot_entity   real(rk)(:,:,:)      rot_c2s of each elset (elset orientations)
  !>   i_line       integer(ik)          line index
  !>   num_params   integer              parameters of the descriptor
  !>   table_name   character(:)         section name, for messages
  subroutine resolve_orientations(mesh, name)

    type(mesh_file_t), intent(inout) :: mesh
    character(len=*), intent(in) :: name

    integer(ik), allocatable :: order(:)
    real(rk), allocatable :: rot_entity(:, :, :)
    integer(ik) :: i_line
    integer :: num_params
    character(len=:), allocatable :: table_name

    allocate (mesh%elt_rot_c2s(3, 3, mesh%num_elts))
    mesh%elt_rot_c2s = 0.0_rk
    mesh%has_ori = mesh%raw_elt_ori%defined .or. mesh%raw_elset_ori%defined
    if (.not. mesh%has_ori) return

    if (mesh%raw_elt_ori%defined) then
      ! Element orientations: ids are 3D element numbers.
      table_name = '$ElementOrientations (' // mesh%raw_elt_ori%source // ')'
      mesh%ori_from_elts = .true.
      call parse_ori_header(mesh%raw_elt_ori%header, table_name, name, mesh%ori_descriptor, &
        & mesh%ori_convention)
      if (allocated(order)) deallocate (order)
      allocate (order(mesh%num_elts))
      order(:) = map_elt_ids(mesh%raw_elt_ori, mesh%num_elts, table_name, name)
      num_params = ori_num_params(mesh%ori_descriptor)
      do i_line = 1_ik, mesh%num_elts
        call check_num_params(mesh%raw_elt_ori, i_line, num_params, table_name, name)
        mesh%elt_rot_c2s(:, :, order(i_line)) = rot_c2s_from_descriptor(mesh%ori_descriptor, &
          & mesh%ori_convention, mesh%raw_elt_ori%vals(1:num_params, i_line))
      end do
    else
      ! Elset orientations: ids are elset ids.
      table_name = '$ElsetOrientations (' // mesh%raw_elset_ori%source // ')'
      call parse_ori_header(mesh%raw_elset_ori%header, table_name, name, mesh%ori_descriptor, &
        & mesh%ori_convention)
      if (allocated(order)) deallocate (order)
      allocate (order(mesh%num_elsets))
      order(:) = map_elset_ids(mesh, mesh%raw_elset_ori, table_name, name)
      num_params = ori_num_params(mesh%ori_descriptor)
      allocate (rot_entity(3, 3, mesh%num_elsets))
      do i_line = 1_ik, mesh%num_elsets
        call check_num_params(mesh%raw_elset_ori, i_line, num_params, table_name, name)
        rot_entity(:, :, order(i_line)) = rot_c2s_from_descriptor(mesh%ori_descriptor, &
          & mesh%ori_convention, mesh%raw_elset_ori%vals(1:num_params, i_line))
      end do
      mesh%elt_rot_c2s = rot_entity(:, :, mesh%elt_elset)
    end if

  end subroutine resolve_orientations

  !> Groups (phases): one per elset; default phase 1.
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_file_t)  whole-mesh content
  !>   name  (in)     character          mesh file name, for messages
  !> Locals:
  !>   order       integer(ik)(:)  elset index of each line
  !>   i_line      integer(ik)     line index
  !>   table_name  character(:)    section name, for messages
  subroutine resolve_groups(mesh, name)

    type(mesh_file_t), intent(inout) :: mesh
    character(len=*), intent(in) :: name

    integer(ik), allocatable :: order(:)
    integer(ik) :: i_line
    character(len=:), allocatable :: table_name

    ! Default: a single phase.
    allocate (mesh%elset_phase(mesh%num_elsets))
    mesh%elset_phase = 1_ik
    if (mesh%raw_groups%defined) then
      table_name = '$Groups (' // mesh%raw_groups%source // ')'
      if (allocated(order)) deallocate (order)
      allocate (order(mesh%num_elsets))
      order(:) = map_elset_ids(mesh, mesh%raw_groups, table_name, name)
      do i_line = 1_ik, mesh%num_elsets
        mesh%elset_phase(order(i_line)) = nint(mesh%raw_groups%vals(1, i_line), ik)
      end do
      if (any(mesh%elset_phase < 1_ik)) call fatal_error('mesh_file_finalize', name // ': ' // &
        & table_name // ': groups (phases) must be positive')
    end if
    mesh%num_phases = maxval(mesh%elset_phase)
    mesh%elt_phase = mesh%elset_phase(mesh%elt_elset)

  end subroutine resolve_groups

  !> Value sections: per elset (by elset index) or per element (by 3D element number).
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_file_t)  whole-mesh content
  !>   name  (in)     character          mesh file name, for messages
  !> Locals:
  !>   i_sec       integer         value-section index
  !>   order       integer(ik)(:)  entity index of each line
  !>   num_lines   integer(ik)     lines of the section
  !>   width       integer(ik)     largest number of values on a line
  !>   i_line      integer(ik)     line index
  !>   table_name  character(:)    section name, for messages
  subroutine resolve_values(mesh, name)

    type(mesh_file_t), intent(inout) :: mesh
    character(len=*), intent(in) :: name

    integer :: i_sec
    integer(ik), allocatable :: order(:)
    integer(ik) :: num_lines
    integer(ik) :: width
    integer(ik) :: i_line
    character(len=:), allocatable :: table_name

    do i_sec = 1, num_value_sections
      associate (raw => mesh%raw_values(i_sec), out => mesh%values(i_sec))
        if (.not. raw%defined) cycle
        table_name = '$' // trim(value_section_names(i_sec)) // ' (' // raw%source // ')'
        if (value_section_is_elset(i_sec)) then
          if (allocated(order)) deallocate (order)
          allocate (order(mesh%num_elsets))
          order(:) = map_elset_ids(mesh, raw, table_name, name)
        else
          if (allocated(order)) deallocate (order)
          allocate (order(mesh%num_elts))
          order(:) = map_elt_ids(raw, mesh%num_elts, table_name, name)
        end if
        num_lines = int(size(raw%ids), ik)
        width = maxval(raw%num)
        out%defined = .true.
        out%source = raw%source
        out%num_values_hdr = parse_header_count(raw%header)
        allocate (out%num_values(num_lines), out%values(width, num_lines))
        out%values = 0.0_rk
        do i_line = 1_ik, num_lines
          out%num_values(order(i_line)) = raw%num(i_line)
          out%values(:, order(i_line)) = raw%vals(1:width, i_line)
        end do
      end associate
    end do

  end subroutine resolve_values

  ! ==========================================================================================
  ! Helpers
  ! ==========================================================================================

  !> Entity index of each line of an elset-id table: every elset exactly once.
  !>
  !> Arguments:
  !>   mesh        (in)  type(mesh_file_t)  whole-mesh content (elsets resolved)
  !>   table       (in)  type(id_table_t)   the table
  !>   table_name  (in)  character          section name, for messages
  !>   name        (in)  character          mesh file name, for messages
  !>   order       (out, result) integer(ik)(:)  elset index of each line
  !> Locals:
  !>   max_id   integer(ik)     largest elset id
  !>   idx_of   integer(ik)(:)  elset index of each id (0: not an elset)
  !>   seen     logical(:)      whether each elset has a line
  !>   i_line   integer(ik)     line index
  function map_elset_ids(mesh, table, table_name, name) result(order)

    type(mesh_file_t), intent(in) :: mesh
    type(id_table_t), intent(in) :: table
    character(len=*), intent(in) :: table_name
    character(len=*), intent(in) :: name
    integer(ik), allocatable :: order(:)

    integer(ik) :: max_id
    integer(ik), allocatable :: idx_of(:)
    logical, allocatable :: seen(:)
    integer(ik) :: i_line

    ! One line per elset.
    if (size(table%ids) /= mesh%num_elsets) call fatal_error('mesh_file_finalize', name // &
      & ': ' // table_name // ' has ' // int_text(int(size(table%ids), ik)) // ' lines for ' // &
      & int_text(mesh%num_elsets) // ' elsets')

    ! Map ids to elset indices; unknown and repeated ids are errors.
    max_id = maxval(mesh%elset_id)
    allocate (idx_of(max_id), seen(mesh%num_elsets), order(size(table%ids)))
    idx_of = 0_ik
    idx_of(mesh%elset_id) = [(i_line, i_line = 1_ik, mesh%num_elsets)]
    seen = .false.
    do i_line = 1_ik, int(size(table%ids), ik)
      associate (id => table%ids(i_line))
        if (id < 1_ik .or. id > max_id) then
          call fatal_error('mesh_file_finalize', name // ': ' // table_name // ': elset ' // &
            & int_text(id) // ' does not exist')
        end if
        if (idx_of(id) == 0_ik) call fatal_error('mesh_file_finalize', name // ': ' // &
          & table_name // ': elset ' // int_text(id) // ' does not exist')
        if (seen(idx_of(id))) call fatal_error('mesh_file_finalize', name // ': ' // &
          & table_name // ': elset ' // int_text(id) // ' given twice')
        seen(idx_of(id)) = .true.
        order(i_line) = idx_of(id)
      end associate
    end do

  end function map_elset_ids

  !> Element number of each line of an element-id table: ids 1..num_elts, each exactly once.
  !>
  !> Arguments:
  !>   table       (in)  type(id_table_t)  the table
  !>   num_elts    (in)  integer(ik)       number of 3D elements
  !>   table_name  (in)  character         section name, for messages
  !>   name        (in)  character         mesh file name, for messages
  !>   order       (out, result) integer(ik)(:)  element number of each line
  !> Locals:
  !>   seen    logical(:)   whether each element has a line
  !>   i_line  integer(ik)  line index
  function map_elt_ids(table, num_elts, table_name, name) result(order)

    type(id_table_t), intent(in) :: table
    integer(ik), intent(in) :: num_elts
    character(len=*), intent(in) :: table_name
    character(len=*), intent(in) :: name
    integer(ik), allocatable :: order(:)

    logical, allocatable :: seen(:)
    integer(ik) :: i_line

    ! One line per element.
    if (size(table%ids) /= num_elts) call fatal_error('mesh_file_finalize', name // ': ' // &
      & table_name // ' has ' // int_text(int(size(table%ids), ik)) // ' lines for ' // &
      & int_text(num_elts) // ' elements')
    allocate (seen(num_elts), order(num_elts))
    seen = .false.
    do i_line = 1_ik, num_elts
      associate (id => table%ids(i_line))
        if (id < 1_ik .or. id > num_elts) call fatal_error('mesh_file_finalize', name // ': ' &
          & // table_name // ': element ' // int_text(id) // ' outside 1..' // &
          & int_text(num_elts))
        if (seen(id)) call fatal_error('mesh_file_finalize', name // ': ' // table_name // &
          & ': element ' // int_text(id) // ' given twice')
        seen(id) = .true.
        order(i_line) = id
      end associate
    end do

  end function map_elt_ids

  !> Descriptor and convention codes of an orientation header token "descriptor:convention".
  !>
  !> Arguments:
  !>   header      (in)   character  the token
  !>   table_name  (in)   character  section name, for messages
  !>   name        (in)   character  mesh file name, for messages
  !>   descriptor  (out)  integer    descriptor code
  !>   convention  (out)  integer    convention code
  !> Locals:
  !>   colon  integer  position of the colon
  subroutine parse_ori_header(header, table_name, name, descriptor, convention)

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: table_name
    character(len=*), intent(in) :: name
    integer, intent(out) :: descriptor
    integer, intent(out) :: convention

    integer :: colon

    ! Split at the colon; both parts must be known names.
    colon = index(header, ':')
    descriptor = 0
    convention = 0
    if (colon > 1) then
      descriptor = ori_descriptor_from_name(header(1:colon - 1))
      convention = ori_convention_from_name(header(colon + 1:))
    end if
    if (descriptor == 0 .or. convention == 0) call fatal_error('mesh_file_finalize', name // &
      & ': ' // table_name // ': invalid orientation descriptor "' // header // '" (expected ' &
      & // 'descriptor:convention, e.g. rodrigues:active)')

  end subroutine parse_ori_header

  !> An orientation line must have exactly the descriptor's number of parameters.
  !>
  !> Arguments:
  !>   table       (in)  type(id_table_t)  the table
  !>   i_line      (in)  integer(ik)       line index
  !>   num_params  (in)  integer           parameters of the descriptor
  !>   table_name  (in)  character         section name, for messages
  !>   name        (in)  character         mesh file name, for messages
  subroutine check_num_params(table, i_line, num_params, table_name, name)

    type(id_table_t), intent(in) :: table
    integer(ik), intent(in) :: i_line
    integer, intent(in) :: num_params
    character(len=*), intent(in) :: table_name
    character(len=*), intent(in) :: name

    if (table%num(i_line) /= int(num_params, ik)) call fatal_error('mesh_file_finalize', name &
      & // ': ' // table_name // ': id ' // int_text(table%ids(i_line)) // ' has ' // &
      & int_text(table%num(i_line)) // ' parameters, expected ' // int_text(int(num_params, ik)))

  end subroutine check_num_params

  !> Read a line holding one non-negative count.
  !>
  !> Arguments:
  !>   file   (inout)  type(text_file_t)  the open file
  !>   what   (in)     character          what is counted, for messages
  !>   count  (out, result) integer(ik)   the count
  !> Locals:
  !>   line    character(:), allocatable  current line
  !>   values  integer(ik)(1)             the count as read
  function read_count(file, what) result(count)

    type(text_file_t), intent(inout) :: file
    character(len=*), intent(in) :: what
    integer(ik) :: count

    character(len=:), allocatable :: line
    integer(ik) :: values(1)

    call text_require_line(file, line, what)
    call parse_ints(file, line, values, what)
    count = values(1)
    if (count < 0_ik) call text_error(file, what // ' must not be negative')

  end function read_count

  !> Read a line holding one label (a single token).
  !>
  !> Arguments:
  !>   file   (inout)  type(text_file_t)  the open file
  !>   what   (in)     character          what the label is, for messages
  !>   label  (out, result) character(:), allocatable
  function read_label(file, what) result(label)

    type(text_file_t), intent(inout) :: file
    character(len=*), intent(in) :: what
    character(len=:), allocatable :: label

    call text_require_line(file, label, what)
    if (count_tokens(label) /= 1) call text_error(file, what // ': expected one word, found "' &
      & // label // '"')
    label = trim(adjustl(label))

  end function read_label

  !> Read the footer line "$End<section>".
  !>
  !> Arguments:
  !>   file     (inout)  type(text_file_t)  the open file
  !>   section  (in)     character          section name (without $)
  !> Locals:
  !>   line  character(:), allocatable  current line
  subroutine read_footer(file, section)

    type(text_file_t), intent(inout) :: file
    character(len=*), intent(in) :: section

    character(len=:), allocatable :: line

    call text_require_line(file, line, '$End' // section)
    if (line /= '$End' // section) call text_error(file, 'expected $End' // section // &
      & ', found "' // line // '" (wrong count or extra lines in $' // section // ')')

  end subroutine read_footer

  !> A mesh-file section may appear only once.
  !>
  !> Arguments:
  !>   file     (in)     type(text_file_t)  the open file
  !>   seen     (inout)  logical            whether the section was read before
  !>   section  (in)     character          section name
  subroutine check_once(file, seen, section)

    type(text_file_t), intent(in) :: file
    logical, intent(inout) :: seen
    character(len=*), intent(in) :: section

    if (seen) call text_error(file, 'section $' // section // ' given twice')
    seen = .true.

  end subroutine check_once

  !> A section read again: an error in the mesh file; from an auxiliary file it replaces the
  !> earlier definition (logged).
  !>
  !> Arguments:
  !>   table    (in)  type(id_table_t)   the section's storage
  !>   file     (in)  type(text_file_t)  the open file
  !>   section  (in)  character          section name
  !>   is_aux   (in)  logical            reading an auxiliary file
  subroutine check_redefinition(table, file, section, is_aux)

    type(id_table_t), intent(in) :: table
    type(text_file_t), intent(in) :: file
    character(len=*), intent(in) :: section
    logical, intent(in) :: is_aux

    if (.not. table%defined) return
    if (is_aux) then
      call log_info('$' // section // ' from `' // file%name // "' replaces the one from `" // &
        & table%source // "'.", 1)
    else
      call text_error(file, 'section $' // section // ' given twice')
    end if

  end subroutine check_redefinition

  !> Fatal error if a global size exceeds huge(ik) (ROADMAP D28).
  !>
  !> Arguments:
  !>   file  (in)  type(text_file_t)  the open file
  !>   size  (in)  integer(i8)        the size
  !>   what  (in)  character          what the size is
  !> Locals:
  !>   text  character(24)  size as text
  subroutine check_size(file, size, what)

    type(text_file_t), intent(in) :: file
    integer(i8), intent(in) :: size
    character(len=*), intent(in) :: what

    character(len=24) :: text

    if (size <= int(huge(1_ik), i8)) return
    write (text, '(i0)') size
    call text_error(file, what // ' ' // trim(text) // ' exceeds the largest integer (' // &
      & int_text(huge(1_ik)) // '): the mesh is too large for the 32-bit indices of FEPX')

  end subroutine check_size

  !> Value-section index of a section name (0 if not a value section).
  !>
  !> Arguments:
  !>   section  (in)  character  section name
  !>   idx      (out, result) integer
  !> Locals:
  !>   i  integer  candidate index
  pure function value_section_index(section) result(idx)

    character(len=*), intent(in) :: section
    integer :: idx

    integer :: i

    idx = 0
    do i = 1, num_value_sections
      if (section == trim(value_section_names(i))) idx = i
    end do

  end function value_section_index

  !> The integer of a value-section header token (checked when read).
  !>
  !> Arguments:
  !>   header  (in)  character  the token
  !>   count   (out, result) integer(ik)
  function parse_header_count(header) result(count)

    character(len=*), intent(in) :: header
    integer(ik) :: count

    read (header, *) count

  end function parse_header_count

  !> Name of a 3D element type, for messages.
  !>
  !> Arguments:
  !>   elt_type  (in)  integer(ik)  Gmsh type code
  !>   name      (out, result) character(:), allocatable
  pure function elt_type_name(elt_type) result(name)

    integer(ik), intent(in) :: elt_type
    character(len=:), allocatable :: name

    select case (elt_type)
    case (elt_type_tet10)
      name = '10-node tetrahedra'
    case (elt_type_hex20)
      name = '20-node hexahedra'
    case default
      name = 'unknown'
    end select

  end function elt_type_name

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

end module mesh_file_mod
