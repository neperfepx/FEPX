! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Dirichlet (velocity) boundary conditions: general conditions on node sets or single nodes,
!> the boundary-condition presets (ROADMAP D14: presets are separate from the control drivers),
!> their resolution to the rank's local dofs, and the gage length.
!>
!> A condition prescribes one velocity component on every node of a node set (or on one node,
!> given by its global id as text). The presets build their conditions for a loading axis a
!> (1 x, 2 y, 3 z) with b, c the next axes in cyclic order (a = z: b = x, c = y); the velocity of
!> the moving face is the strain rate times the gage length (the extent of the mesh along a):
!>   uniaxial_grip      face a0: all three components 0; face a1: b, c 0, a moving
!>   uniaxial_symmetry  face b0: b 0; face c0: c 0; face a0: a 0; face a1: a moving
!>   uniaxial_minimal   face a0: a 0; face a1: a moving; corner x0y0z0: b 0 and c 0; corner with
!>                      coordinate 1 along b (0 along a and c): c 0
!>   triaxial           faces x0, x1: x; y0, y1: y; z0, z1: z (values set by the control driver,
!>                      zero here); with velocities: bc_triaxial_specs (step 5.3, D45)
!> Corners are Neper's node sets "x<i>y<j>z<k>".
!>
!> Source : port of v2 boundary_conditions_mod2_legacy.f90 (calc_bcs_uniaxial_grip,
!>          calc_bcs_uniaxial_symmetry, calc_bcs_uniaxial_minimal with all five conditions,
!>          CHANGES.md C25, calc_bcs_triaxial; oriented along the loading direction as v2-fixed
!>          does after C26, the same conditions as v2's strain-rate "minimal" after C27) and
!>          boundary_conditions_mod3.f90::calc_bcs_general_vel (node set or single node, conflict
!>          check); gage length: v2 mesh_type_mod.f90::mesh_lengths (bounding-box extent).
!>
!> Public parameters:
!>   bc_preset_grip, bc_preset_symmetry, bc_preset_minimal, bc_preset_triaxial  preset codes
!>
!> Public types:
!>   bc_spec_t    one condition (components documented at the type)
!>
!> Public routines:
!>   bc_preset_specs   the conditions of a preset
!>   bc_triaxial_specs the conditions of the triaxial preset under strain control (D45)
!>   bc_resolve        prescribed-dof mask and values on the rank's local dofs
!>   bc_gage_length    extent of the mesh along an axis (all ranks)
!>   bc_lengths        extents of given coordinates along the three axes (all ranks; step 5.2)
module bc_mod

  use kinds_mod, only: rk, ik
  use error_mod, only: fatal_error, internal_error
  use mesh_mod, only: mesh_t, node_loc_from_glob
  use parallel_mod, only: par_min, par_max

  implicit none

  private

  public :: bc_preset_grip, bc_preset_symmetry, bc_preset_minimal, bc_preset_triaxial
  public :: bc_spec_t, bc_preset_specs, bc_triaxial_specs, bc_resolve, bc_gage_length, &
    & bc_lengths

  ! Preset codes.
  integer, parameter :: bc_preset_grip = 1
  integer, parameter :: bc_preset_symmetry = 2
  integer, parameter :: bc_preset_minimal = 3
  integer, parameter :: bc_preset_triaxial = 4

  !> One velocity condition.
  !>
  !> Components:
  !>   nset   character(:), allocatable  node-set label (e.g. "z1", "x0y0z0"), or a global node
  !>                                     id written as an integer
  !>   comp   integer                    velocity component (1 x, 2 y, 3 z)
  !>   value  real(rk)                   prescribed velocity [mm/s]
  type :: bc_spec_t
    character(len=:), allocatable :: nset
    integer :: comp = 0
    real(rk) :: value = 0.0_rk
  end type bc_spec_t

  ! Axis letters.
  character(len=1), parameter :: axis_letter(3) = ['x', 'y', 'z']

contains

  !> The conditions of a preset for loading along axis a.
  !>
  !> Arguments:
  !>   preset    (in)  integer           preset code (bc_preset_*)
  !>   axis      (in)  integer           loading axis a (1 x, 2 y, 3 z)
  !>   vel_load  (in)  real(rk)          velocity of the moving face (strain rate times gage
  !>                                     length) [mm/s]
  !>   specs     (out) type(bc_spec_t)(:), allocatable  conditions
  !> Locals:
  !>   b, c      integer                 the next axes after a in cyclic order
  !>   corner_b  character(6)            corner with coordinate 1 along b
  subroutine bc_preset_specs(preset, axis, vel_load, specs)

    integer, intent(in) :: preset
    integer, intent(in) :: axis
    real(rk), intent(in) :: vel_load
    type(bc_spec_t), allocatable, intent(out) :: specs(:)

    integer :: b, c
    character(len=6) :: corner_b

    ! Axes and the corner of uniaxial_minimal.
    if (axis < 1 .or. axis > 3) call internal_error('bc_preset_specs', 'the loading axis must be ' &
      & // '1, 2 or 3 (x, y or z)')
    b = mod(axis, 3) + 1
    c = mod(axis + 1, 3) + 1
    corner_b = 'x0y0z0'
    corner_b(2*b:2*b) = '1'

    ! Conditions of the preset (v2's order, oriented along a). Assigned one by one: array
    ! constructors of this type (allocatable character component) lost their contents with
    ! gfortran 13 -O3 (ROADMAP section G, step 2.2).
    select case (preset)
    case (bc_preset_grip)
      allocate (specs(6))
      call set_spec(specs(1), face(axis, 0), b, 0.0_rk)
      call set_spec(specs(2), face(axis, 0), c, 0.0_rk)
      call set_spec(specs(3), face(axis, 0), axis, 0.0_rk)
      call set_spec(specs(4), face(axis, 1), b, 0.0_rk)
      call set_spec(specs(5), face(axis, 1), c, 0.0_rk)
      call set_spec(specs(6), face(axis, 1), axis, vel_load)
    case (bc_preset_symmetry)
      allocate (specs(4))
      call set_spec(specs(1), face(b, 0), b, 0.0_rk)
      call set_spec(specs(2), face(c, 0), c, 0.0_rk)
      call set_spec(specs(3), face(axis, 0), axis, 0.0_rk)
      call set_spec(specs(4), face(axis, 1), axis, vel_load)
    case (bc_preset_minimal)
      allocate (specs(5))
      call set_spec(specs(1), face(axis, 0), axis, 0.0_rk)
      call set_spec(specs(2), face(axis, 1), axis, vel_load)
      call set_spec(specs(3), 'x0y0z0', b, 0.0_rk)
      call set_spec(specs(4), 'x0y0z0', c, 0.0_rk)
      call set_spec(specs(5), corner_b, c, 0.0_rk)
    case (bc_preset_triaxial)
      allocate (specs(6))
      call set_spec(specs(1), 'x0', 1, 0.0_rk)
      call set_spec(specs(2), 'x1', 1, 0.0_rk)
      call set_spec(specs(3), 'y0', 2, 0.0_rk)
      call set_spec(specs(4), 'y1', 2, 0.0_rk)
      call set_spec(specs(5), 'z0', 3, 0.0_rk)
      call set_spec(specs(6), 'z1', 3, 0.0_rk)
    case default
      call internal_error('bc_preset_specs', 'unknown boundary-condition preset')
    end select

  end subroutine bc_preset_specs

  !> Conditions of the triaxial preset under strain control (D45): faces x0, y0, z0 with zero
  !> normal velocity (symmetry); each prescribed face x1, y1, z1 with its normal velocity
  !> (strain rate times the original length along the axis); the face x1, y1 or z1 of a free
  !> component has no condition (zero traction).
  !>
  !> Source : port of v2 boundary_conditions_mod2_legacy.f90::calc_bcs_triaxial (normal
  !>          conditions on the six faces, x0 y0 z0 fixed and x1 y1 z1 moving, as v2's triaxial
  !>          drivers set them), with free faces for components without rows (D45).
  !>
  !> Arguments:
  !>   vel         (in)  real(rk)(3)          normal velocity of faces x1, y1, z1 [mm/s]
  !>   prescribed  (in)  logical(3)           faces x1, y1, z1 with a condition
  !>   specs       (out) type(bc_spec_t)(:)   the conditions
  !> Locals:
  !>   axis, n     integer                    axis, number of conditions set
  subroutine bc_triaxial_specs(vel, prescribed, specs)

    real(rk), intent(in) :: vel(3)
    logical, intent(in) :: prescribed(3)
    type(bc_spec_t), allocatable, intent(out) :: specs(:)

    integer :: axis, n

    ! Symmetry faces, then the prescribed moving faces (assigned one by one: see
    ! bc_preset_specs).
    allocate (specs(3 + count(prescribed)))
    n = 0
    do axis = 1, 3
      n = n + 1
      call set_spec(specs(n), face(axis, 0), axis, 0.0_rk)
    end do
    do axis = 1, 3
      if (.not. prescribed(axis)) cycle
      n = n + 1
      call set_spec(specs(n), face(axis, 1), axis, vel(axis))
    end do

  end subroutine bc_triaxial_specs

  !> Prescribed-dof mask and values on the rank's local dofs (owned and ghost), from a list of
  !> conditions. A node set unknown to the mesh, a component outside 1..3, or one dof given two
  !> different values is a fatal error (v2 calc_bcs_general_vel).
  !>
  !> Arguments:
  !>   mesh     (in)  type(mesh_t)       rank-local mesh (nsets with local node numbers)
  !>   specs    (in)  type(bc_spec_t)(:) conditions
  !>   bc_mask  (out) logical(:)         (3 num_nodes_loc) prescribed dofs
  !>   bc_vel   (out) real(rk)(:)        (3 num_nodes_loc) prescribed values (0 elsewhere) [mm/s]
  !> Locals:
  !>   i_spec, i_set, i  integer          loop indices
  !>   nodes    integer(ik)(:)            local nodes of the condition
  !>   dof      integer(ik)               local dof
  !>   id_glob, ios      integer          global node id of a numeric label, read status
  !>   found    logical                   label found
  subroutine bc_resolve(mesh, specs, bc_mask, bc_vel)

    type(mesh_t), intent(in) :: mesh
    type(bc_spec_t), intent(in) :: specs(:)
    logical, intent(out) :: bc_mask(:)
    real(rk), intent(out) :: bc_vel(:)

    integer :: i_spec, i_set, i, ios
    integer(ik), allocatable :: nodes(:)
    integer(ik) :: dof, id_glob
    logical :: found

    bc_mask = .false.
    bc_vel = 0.0_rk
    do i_spec = 1, size(specs)
      ! Component.
      if (specs(i_spec)%comp < 1 .or. specs(i_spec)%comp > 3) then
        call internal_error('bc_resolve', 'velocity component of node set ' // specs(i_spec)%nset &
          & // ' must be 1, 2 or 3 (x, y or z)')
      end if

      ! Nodes: a global node id, or a node set.
      read (specs(i_spec)%nset, *, iostat=ios) id_glob
      if (ios == 0) then
        nodes = [node_loc_from_glob(mesh, id_glob)]
        nodes = pack(nodes, nodes > 0_ik)
      else
        found = .false.
        do i_set = 1, size(mesh%nsets)
          if (mesh%nsets(i_set)%label == specs(i_spec)%nset) then
            nodes = mesh%nsets(i_set)%nodes
            found = .true.
            exit
          end if
        end do
        if (.not. found) call fatal_error('bc_resolve', 'node set ' // specs(i_spec)%nset // &
          & ' is not in the mesh')
      end if

      ! Prescribe the component on every local node; a dof given two values is an error.
      do i = 1, size(nodes)
        dof = 3_ik*(nodes(i) - 1_ik) + int(specs(i_spec)%comp, ik)
        if (bc_mask(dof) .and. abs(bc_vel(dof) - specs(i_spec)%value) > 0.0_rk) then
          call fatal_error('bc_resolve', 'conflicting velocities prescribed on a node of node ' &
            & // 'set ' // specs(i_spec)%nset // ' (component ' // axis_letter( &
            & specs(i_spec)%comp) // ')')
        end if
        bc_mask(dof) = .true.
        bc_vel(dof) = specs(i_spec)%value
      end do
    end do

  end subroutine bc_resolve

  !> Extent of the mesh along an axis (bounding box over all ranks), v2's gage length.
  !>
  !> Arguments:
  !>   mesh    (in)  type(mesh_t)  rank-local mesh
  !>   axis    (in)  integer       axis (1 x, 2 y, 3 z)
  !>   length  (out, result) real(rk)  max - min of the coordinate over all nodes [mm]
  !> Locals:
  !>   lengths real(rk)(3)         extents along the three axes [mm]
  function bc_gage_length(mesh, axis) result(length)

    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: axis
    real(rk) :: length

    real(rk) :: lengths(3)

    ! Extents of the original coordinates.
    lengths = bc_lengths(mesh, mesh%coo_s)
    length = lengths(axis)

  end function bc_gage_length

  !> Extents of coordinates along the three axes (bounding box over all ranks): the
  !> face-to-face lengths of a box-shaped domain, for the macroscopic strains of the run limits
  !> (D44). Collective.
  !>
  !> Source : port of v2 mesh_type_mod.f90::mesh_lengths (bounding box of the current
  !>          coordinates).
  !>
  !> Arguments:
  !>   mesh     (in)  type(mesh_t)    rank-local mesh
  !>   coo_s    (in)  real(rk)(:, :)  (3, num_nodes_loc) coordinates of the local nodes [mm]
  !>   lengths  (out, result) real(rk)(3)  max - min of each coordinate over all nodes [mm]
  !> Locals:
  !>   lo, hi          real(rk)(3)  local minima and maxima [mm]
  !>   lo_all, hi_all  real(rk)     minimum and maximum of one axis over the ranks [mm]
  !>   axis            integer      axis
  function bc_lengths(mesh, coo_s) result(lengths)

    type(mesh_t), intent(in) :: mesh
    real(rk), intent(in) :: coo_s(:, :)
    real(rk) :: lengths(3)

    real(rk) :: lo(3), hi(3), lo_all, hi_all
    integer :: axis

    ! Local extents (a rank without nodes contributes nothing), then over the ranks.
    lo = huge(1.0_rk)
    hi = -huge(1.0_rk)
    if (mesh%num_nodes_loc > 0_ik) then
      lo = minval(coo_s(:, 1:mesh%num_nodes_loc), dim=2)
      hi = maxval(coo_s(:, 1:mesh%num_nodes_loc), dim=2)
    end if
    do axis = 1, 3
      call par_max(hi(axis), hi_all)
      call par_min(lo(axis), lo_all)
      lengths(axis) = hi_all - lo_all
    end do

  end function bc_lengths

  !> Name of face side (0 or 1) of an axis, e.g. face(3, 1) = "z1".
  !>
  !> Arguments:
  !>   axis  (in)  integer  axis (1, 2, 3)
  !>   side  (in)  integer  0 or 1
  pure function face(axis, side) result(name)

    integer, intent(in) :: axis
    integer, intent(in) :: side
    character(len=2) :: name

    name = axis_letter(axis) // achar(iachar('0') + side)

  end function face

  !> Set one condition.
  !>
  !> Arguments:
  !>   bc     (inout) type(bc_spec_t)  condition
  !>   nset   (in)    character        node-set label
  !>   comp   (in)    integer          component
  !>   value  (in)    real(rk)         velocity [mm/s]
  pure subroutine set_spec(bc, nset, comp, value)

    type(bc_spec_t), intent(inout) :: bc
    character(len=*), intent(in) :: nset
    integer, intent(in) :: comp
    real(rk), intent(in) :: value

    bc%nset = nset
    bc%comp = comp
    bc%value = value

  end subroutine set_spec

end module bc_mod
