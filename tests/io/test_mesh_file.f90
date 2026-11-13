! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of mesh_file_mod (Neper mesh parser) on the fixtures of tests/io/fixtures and on two
!> Neper-written regression meshes.
!>
!> Source : FEPX (ROADMAP step 2.2). Expected values are those written in the fixture files
!>          (tests/io/fixtures/README.md) and in the Neper meshes.
!>
!> Arguments (command line): 1. fixture directory (tests/io/fixtures); 2. regression mesh
!> directory (tests/regression/cases/meshes).
!>
!> Checks:
!>   - tiny.msh: every section (nodes, elements with skipped 0D/2D elements and nodes in the
!>     file's Gmsh order, elset numbering by increasing id, partitions, nsets, fasets mapped to
!>     3D element numbers with each face line kept in Gmsh order (D107), periodicity,
!>     elset orientations out of order converted with orientation_mod, groups), the warning for the
!>     unknown section and the summary line of the log;
!>   - tiny.msh + tiny.opt + tiny.ori (+ tiny.opt again): value sections by elset index / element
!>     number (out of order, variable value counts), element orientations used instead of the
!>     elset orientations, the replacement of a section from an auxiliary file (logged);
!>   - poly_part4.msh (Neper -part 4): sizes, first element, elsets, and the partition-tag layout
!>     (ROADMAP 2.2): partition = third tag, 114/113/114/113 elements, node partitions
!>     198/198/198/198, elements and nodes not contiguous by partition; every faset face lies on
!>     its element (checked by the parser);
!>   - sx2_part4.msh: every node except 248 in a partition of one of its elements (that one is
!>     in partition 2, its elements elsewhere);
!>   - sx2.msh: unpartitioned (all tags 0, no $NodePartitions); 2 elsets with the same
!>     orientation (a single crystal meshed as two cells).
!>
!> Locals:
!>   fixtures   character(512)  fixture directory
!>   meshes     character(512)  regression mesh directory
!>   mesh       mesh_file_t     parsed mesh
!>   unit_log   integer         scratch unit receiving the log
!>   log_text   character(:)    the log, concatenated
!>   expected   real(rk)(3,3)   expected rotation
!>   i_part     integer(ik)     partition index
!>   i_node     integer(ik)     node index
!>   i_elt      integer(ik)     element index
!>   ok         logical         result of a loop check
!>   aux        character(512)(3)  auxiliary file names
program test_mesh_file

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk, ik
  use log_mod, only: log_set_unit
  use orientation_mod, only: rot_c2s_from_descriptor, ori_des_euler_bunge, ori_des_quaternion, &
    & ori_des_rodrigues, ori_conv_active, ori_conv_passive
  use mesh_file_mod, only: mesh_file_t, mesh_file_load, elt_type_tet10
  use assert_mod, only: assert_equal, assert_close, assert_true, test_summary

  implicit none

  character(len=512) :: fixtures
  character(len=512) :: meshes
  type(mesh_file_t) :: mesh
  integer :: unit_log
  character(len=:), allocatable :: log_text
  real(rk) :: expected(3, 3)
  integer(ik) :: i_part
  integer(ik) :: i_node
  integer(ik) :: i_elt
  logical :: ok
  character(len=512) :: aux(3)

  call get_command_argument(1, fixtures)
  call get_command_argument(2, meshes)

  ! --- tiny.msh -------------------------------------------------------------------------------

  call capture_log_start(unit_log)
  call mesh_file_load(mesh, trim(fixtures) // '/tiny.msh')
  log_text = capture_log_end(unit_log)

  ! Header sections.
  call assert_equal(mesh%version, '2.4', 'tiny version')
  call assert_equal(mesh%domain, 'cube', 'tiny domain')

  ! Nodes.
  call assert_equal(mesh%num_nodes, 14_ik, 'tiny num_nodes')
  call assert_close(mesh%coo(:, 5), [1.0_rk, 1.0_rk, 1.0_rk], 'tiny coo 5')
  call assert_close(mesh%coo(:, 14), [0.5_rk, 0.5_rk, 1.0_rk], 'tiny coo 14')

  ! Elements: the 0D and 2D elements are skipped; node order kept (Gmsh order, ROADMAP D35).
  call assert_equal(mesh%elt_type, elt_type_tet10, 'tiny elt_type')
  call assert_equal(mesh%num_elt_nodes, 10_ik, 'tiny num_elt_nodes')
  call assert_equal(mesh%num_elts, 2_ik, 'tiny num_elts')
  call assert_equal(mesh%elt_id_first, 3_ik, 'tiny elt_id_first')
  call assert_int_array(mesh%elt_nodes(:, 1), [1, 2, 3, 4, 6, 7, 8, 9, 10, 11], 'tiny elt 1')
  call assert_int_array(mesh%elt_nodes(:, 2), [2, 3, 4, 5, 7, 10, 11, 12, 14, 13], 'tiny elt 2')

  ! Elsets 3 and 7 numbered 1 and 2.
  call assert_equal(mesh%num_elsets, 2_ik, 'tiny num_elsets')
  call assert_int_array(mesh%elset_id, [3, 7], 'tiny elset_id')
  call assert_int_array(mesh%elt_elset, [1, 2], 'tiny elt_elset')

  ! Partitions.
  call assert_equal(mesh%num_parts, 2_ik, 'tiny num_parts')
  call assert_int_array(mesh%elt_part, [1, 2], 'tiny elt_part')
  call assert_int_array(mesh%node_part, [1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2], &
    & 'tiny node_part')

  ! Node sets.
  call assert_equal(size(mesh%nsets), 2, 'tiny num nsets')
  call assert_equal(mesh%nsets(1)%label, 'face', 'tiny nset 1 label')
  call assert_int_array(mesh%nsets(1)%nodes, [2, 3, 4], 'tiny nset 1 nodes')
  call assert_equal(mesh%nsets(2)%label, 'apex', 'tiny nset 2 label')
  call assert_int_array(mesh%nsets(2)%nodes, [5], 'tiny nset 2 nodes')

  ! Face set: file element id 3 is 3D element 1.
  call assert_equal(size(mesh%fasets), 1, 'tiny num fasets')
  call assert_equal(mesh%fasets(1)%label, 'z0', 'tiny faset label')
  call assert_equal(mesh%fasets(1)%num_face_nodes, 6_ik, 'tiny faset nodes per face')
  call assert_int_array(mesh%fasets(1)%elts, [1], 'tiny faset elts')
  ! The Gmsh-order triangle 1 3 2 of the file (normal -z, out of the body at z = 0) with edge
  ! nodes 8 (1-3), 7 (3-2), 6 (2-1).
  call assert_int_array(mesh%fasets(1)%nodes(:, 1), [1, 3, 2, 8, 7, 6], 'tiny faset face')

  ! Periodicity.
  call assert_equal(mesh%num_periodicity, 1_ik, 'tiny num_periodicity')
  call assert_int_array(mesh%per_secondary, [5], 'tiny per secondary')
  call assert_int_array(mesh%per_primary, [1], 'tiny per primary')
  call assert_int_array(mesh%per_vect(:, 1), [1, 1, 1], 'tiny per vect')

  ! Elset orientations (lines out of order): element 1 is elset 3, element 2 elset 7.
  call assert_true(mesh%has_ori, 'tiny has_ori')
  call assert_true(.not. mesh%ori_from_elts, 'tiny ori from elsets')
  call assert_equal(mesh%ori_descriptor, ori_des_euler_bunge, 'tiny ori descriptor')
  call assert_equal(mesh%ori_convention, ori_conv_active, 'tiny ori convention')
  expected = rot_c2s_from_descriptor(ori_des_euler_bunge, ori_conv_active, &
    & [40.0_rk, 50.0_rk, 60.0_rk])
  call assert_close(mesh%elt_rot_c2s(:, :, 1), expected, 'tiny rot_c2s 1', abs_tol=1e-15_rk)
  expected = rot_c2s_from_descriptor(ori_des_euler_bunge, ori_conv_active, &
    & [10.0_rk, 20.0_rk, 30.0_rk])
  call assert_close(mesh%elt_rot_c2s(:, :, 2), expected, 'tiny rot_c2s 2', abs_tol=1e-15_rk)

  ! Groups.
  call assert_equal(mesh%num_phases, 2_ik, 'tiny num_phases')
  call assert_int_array(mesh%elset_phase, [1, 2], 'tiny elset_phase')
  call assert_int_array(mesh%elt_phase, [1, 2], 'tiny elt_phase')

  ! No value sections in the mesh file.
  call assert_true(.not. any(mesh%values(:)%defined), 'tiny no value sections')

  ! Log: unknown section skipped with a warning; the file lines (step 9.7).
  call assert_true(index(log_text, 'Warning:   - Skipping unknown mesh section ' // &
    & '$FutureSection of `' // trim(fixtures) // "/tiny.msh'.") > 0, &
    & 'tiny warning unknown section')
  call assert_true(index(log_text, 'Info   :   [i] Parsing file `' // trim(fixtures) // &
    & "/tiny.msh'..." // new_line('a') // 'Warning:') > 0, 'tiny parsing line')
  call assert_true(index(log_text, 'Info   :   [i] Parsed file `' // trim(fixtures) // &
    & "/tiny.msh'.") > 0, 'tiny parsed line')

  ! --- tiny.msh + auxiliary files --------------------------------------------------------------

  aux(1) = trim(fixtures) // '/tiny.opt'
  aux(2) = trim(fixtures) // '/tiny.ori'
  aux(3) = trim(fixtures) // '/tiny.opt'
  call capture_log_start(unit_log)
  call mesh_file_load(mesh, trim(fixtures) // '/tiny.msh', aux)
  log_text = capture_log_end(unit_log)

  ! $ElsetCrss: elset 3 (index 1) has 2 values, elset 7 (index 2) one.
  call assert_true(mesh%values(1)%defined, 'aux ElsetCrss defined')
  call assert_equal(mesh%values(1)%num_values_hdr, 12_ik, 'aux ElsetCrss header')
  call assert_int_array(mesh%values(1)%num_values, [2, 1], 'aux ElsetCrss counts')
  call assert_close(mesh%values(1)%values(:, 1), [100.0_rk, 110.0_rk], 'aux ElsetCrss elset 3')
  call assert_close(mesh%values(1)%values(1, 2), 200.0_rk, 'aux ElsetCrss elset 7')

  ! $ElementCrssSat: by element number.
  call assert_true(mesh%values(5)%defined, 'aux ElementCrssSat defined')
  call assert_close(mesh%values(5)%values(1, :), [330.0_rk, 335.0_rk], 'aux ElementCrssSat')
  call assert_true(.not. mesh%values(2)%defined, 'aux ElsetCrssSat undefined')

  ! $ElsetCrySym of tiny.ori accepted: the mesh file's (step 8.5).
  call assert_true(mesh%elset_crysym == 'triclinic', 'aux ElsetCrySym kept')

  ! Element orientations are used instead of the elset orientations.
  call assert_true(mesh%ori_from_elts, 'aux ori from elements')
  call assert_equal(mesh%ori_descriptor, ori_des_quaternion, 'aux ori descriptor')
  call assert_equal(mesh%ori_convention, ori_conv_passive, 'aux ori convention')
  expected = rot_c2s_from_descriptor(ori_des_quaternion, ori_conv_passive, &
    & [1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk])
  call assert_close(mesh%elt_rot_c2s(:, :, 1), expected, 'aux rot_c2s 1', abs_tol=1e-15_rk)
  expected = rot_c2s_from_descriptor(ori_des_quaternion, ori_conv_passive, &
    & [0.5_rk, 0.5_rk, 0.5_rk, 0.5_rk])
  call assert_close(mesh%elt_rot_c2s(:, :, 2), expected, 'aux rot_c2s 2', abs_tol=1e-15_rk)

  ! The second tiny.opt replaced the first (logged).
  call assert_true(index(log_text, '$ElsetCrss from `' // trim(aux(3)) // "' replaces the " // &
    & 'one from `' // trim(aux(1)) // "'.") > 0, 'aux replacement logged')

  ! --- poly_part4.msh: the partition-tag layout of a Neper -part 4 mesh ------------------------

  call capture_log_start(unit_log)
  call mesh_file_load(mesh, trim(meshes) // '/poly_part4.msh')
  log_text = capture_log_end(unit_log)

  call assert_equal(mesh%num_nodes, 792_ik, 'poly4 num_nodes')
  call assert_equal(mesh%num_elts, 454_ik, 'poly4 num_elts')
  call assert_equal(mesh%elt_id_first, 490_ik, 'poly4 elt_id_first')
  call assert_int_array(mesh%elt_nodes(:, 1), [45, 51, 55, 93, 129, 130, 131, 132, 133, 134], &
    & 'poly4 elt 1')
  call assert_equal(mesh%num_elsets, 8_ik, 'poly4 num_elsets')
  call assert_int_array(mesh%elset_id, [1, 2, 3, 4, 5, 6, 7, 8], 'poly4 elset_id')
  call assert_equal(size(mesh%nsets), 26, 'poly4 nsets')
  call assert_equal(size(mesh%fasets), 6, 'poly4 fasets')
  call assert_equal(mesh%ori_descriptor, ori_des_rodrigues, 'poly4 descriptor')
  call assert_equal(mesh%ori_convention, ori_conv_passive, 'poly4 convention')

  ! Partition = third tag; 114, 113, 114 and 113 elements; node partitions as in the file.
  call assert_equal(mesh%num_parts, 4_ik, 'poly4 num_parts')
  call assert_int_array([(int(count(mesh%elt_part == i_part)), i_part = 1_ik, 4_ik)], &
    & [114, 113, 114, 113], 'poly4 elts per part')
  call assert_int_array([(int(count(mesh%node_part == i_part)), i_part = 1_ik, 4_ik)], &
    & [198, 198, 198, 198], 'poly4 nodes per part')

  ! Elements and nodes are not numbered by partition (Neper no longer renumbers them, D108).
  call assert_true(.not. all(mesh%elt_part(2:) >= mesh%elt_part(:mesh%num_elts - 1)), &
    & 'poly4 elements not contiguous by partition')
  call assert_true(.not. all(mesh%node_part(2:) >= mesh%node_part(:mesh%num_nodes - 1)), &
    & 'poly4 nodes not contiguous by partition')

  ! Phases: no $Groups, a single phase.
  call assert_equal(mesh%num_phases, 1_ik, 'poly4 single phase')

  ! --- sx2_part4.msh: a node in a partition without any of its elements -----------------------

  ! Neper's partitioner does not always give a node to a partition that contains one of its
  ! elements: here exactly node 248 is in partition 2, while its elements are all in other
  ! partitions (the distribution must allow an owned node without local elements). (Until step
  ! 9.20, node 945 of poly_part4.msh; its Neper 5.0.1-731-devel mesh has no such node.)
  call capture_log_start(unit_log)
  call mesh_file_load(mesh, trim(meshes) // '/sx2_part4.msh')
  log_text = capture_log_end(unit_log)
  call assert_equal(mesh%num_parts, 4_ik, 'sx2_4 num_parts')
  ok = .true.
  do i_node = 1_ik, mesh%num_nodes
    if (i_node == 248_ik) then
      if (mesh%node_part(i_node) /= 2_ik) ok = .false.
      if (any(mesh%elt_part == 2_ik .and. any(mesh%elt_nodes == i_node, dim=1))) ok = .false.
    else
      if (.not. any(mesh%elt_part == mesh%node_part(i_node) .and. &
        & any(mesh%elt_nodes == i_node, dim=1))) ok = .false.
    end if
  end do
  call assert_true(ok, 'sx2_4 node partitions: node 248 has no element in its own')

  ! --- sx2.msh: unpartitioned, a single crystal meshed as two cells ----------------------------

  call capture_log_start(unit_log)
  call mesh_file_load(mesh, trim(meshes) // '/sx2.msh')
  log_text = capture_log_end(unit_log)
  call assert_equal(mesh%num_elts, 145_ik, 'sx2 num_elts')
  call assert_equal(mesh%num_elsets, 2_ik, 'sx2 num_elsets')
  call assert_equal(mesh%num_parts, 0_ik, 'sx2 unpartitioned')
  call assert_true(all(mesh%node_part == 0_ik), 'sx2 node_part 0')
  call assert_true(all(mesh%elt_part == 0_ik), 'sx2 elt_part 0')
  ok = .true.
  do i_elt = 1_ik, mesh%num_elts
    if (any(abs(mesh%elt_rot_c2s(:, :, i_elt) - mesh%elt_rot_c2s(:, :, 1)) > 0.0_rk)) ok = .false.
  end do
  call assert_true(ok, 'sx2 one orientation for every element')

  call test_summary()

contains

  !> Assert that an integer(ik) array equals the given default-integer values.
  !>
  !> Arguments:
  !>   actual    (in)  integer(ik)(:)  array under test
  !>   expected  (in)  integer(:)      expected values
  !>   label     (in)  character       check name
  !> Locals:
  !>   i  integer  index
  subroutine assert_int_array(actual, expected, label)

    integer(ik), intent(in) :: actual(:)
    integer, intent(in) :: expected(:)
    character(len=*), intent(in) :: label

    integer :: i

    call assert_equal(size(actual), size(expected), label // ' (size)')
    if (size(actual) /= size(expected)) return
    do i = 1, size(actual)
      call assert_equal(int(actual(i)), expected(i), label)
    end do

  end subroutine assert_int_array

  !> Start capturing the log in a scratch file.
  !>
  !> Arguments:
  !>   unit  (out)  integer  scratch unit
  subroutine capture_log_start(unit)

    integer, intent(out) :: unit

    open (newunit=unit, status='scratch', action='readwrite', form='formatted')
    call log_set_unit(unit)

  end subroutine capture_log_start

  !> Stop capturing the log and return it (lines joined with newlines); echo it to standard
  !> output for the test log.
  !>
  !> Arguments:
  !>   unit  (in)  integer  scratch unit
  !>   text  (out, result) character(:), allocatable
  !> Locals:
  !>   line  character(1024)  one line
  !>   ios   integer          iostat of the reads
  function capture_log_end(unit) result(text)

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
      write (output_unit, '(a)') trim(line)
      ! A continuation line of a wrapped message (log_mod, step 9.7: the prefix, then two or
      ! more blanks and no depth or file marker) joins its message with one blank.
      if (len(text) > 0 .and. line(10:11) == '  ') then
        if (verify(line(10:), ' ') > 0) then
          if (index('->.[', line(9 + verify(line(10:), ' '):9 + verify(line(10:), ' '))) == 0) &
            & then
            text = text(1:len(text) - 1) // ' ' // trim(adjustl(line(10:))) // new_line('a')
            cycle
          end if
        end if
      end if
      text = text // trim(line) // new_line('a')
    end do
    close (unit)

  end function capture_log_end

end program test_mesh_file
