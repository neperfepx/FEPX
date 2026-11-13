# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# FEPX test suite.
#
# Layout:
#   framework/  in-house assertion module (library fepx_testing) and its self-tests
#   core/       kinds, constants, error handling and logging (src/core/)
#   math/       tensor algebra (src/math/), property tests against explicit 3x3 algebra
#   material/   crystal structures, slip systems, stiffness (src/material/)
#   element/    element library: shape functions, quadrature (src/element/)
#   fe/         element matrices, EBE operations, Jacobi-PCG, Dirichlet conditions (src/fe/)
#   io/         text files and the Neper mesh parser (src/io/)
#   parallel/   MPI wrappers (src/parallel/), run on 1, 2 and 4 ranks
#   dev/tools/      unit tests of the development tools in dev/tools/ (Python)
#   regression/ regression cases against FEPX's own reference data
#   verification/ cases checked against analytic solutions
#
# Helpers fepx_add_unit_test / fepx_add_mpi_test: src/cmake/FEPXTesting.cmake.

# --- Time limits of the tests that run simulations -----------------------------------------

# The limits given in the test definitions (60, 120, 300, 600 and 3600 s; 3600: the v2 tutorial
# cases, step 8.5) are for an optimized
# build. A Debug build (no optimization, run-time checks) is 2-3 times slower on the finite
# element cases (measured 2026-10-16, step 7.3, Release in dev/tools/ci_local.sh against Debug
# alone on the machine: the sub-incrementation app run 42 s / 88 s, km_poly_fcc_rho_opt
# 66 s / 149 s, irr_sx_bcc 148 s / 262 s, irr_poly_bcc_opt 175 s / 311 s), and the
# parallel test runs of dev/tools/ci_local.sh (ctest -j) slow each case down further, so in Debug
# every limit is multiplied by fepx_timeout_scale. The Checked build (-O2 and Debug's run-time
# checks) is 1.2-1.6 times slower than Release (measured 2026-10-24, step 9.6, ROADMAP D88): its
# limits are doubled. The variables fepx_timeout_<seconds> hold the scaled limits used by the
# test definitions here and in verification/ and regression/.
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(fepx_timeout_scale 4)
elseif(CMAKE_BUILD_TYPE STREQUAL "Checked")
  set(fepx_timeout_scale 2)
else()
  set(fepx_timeout_scale 1)
endif()
foreach(limit IN ITEMS 60 120 300 600 3600)
  math(EXPR fepx_timeout_${limit} "${limit} * ${fepx_timeout_scale}")
endforeach()

# --- Test tiers (step 9.6, ROADMAP D88) ------------------------------------------------------

# Labels "quick" (the unit tests and the tests named in quick_tests.txt), "full" (every test but
# the tutorial cases) and "tutorial"; added at the end of each test directory
# (fepx_apply_test_tiers, src/cmake/FEPXTesting.cmake).
fepx_read_quick_tests(${CMAKE_CURRENT_SOURCE_DIR}/quick_tests.txt)
# At the end of this directory (after verification/ and regression/, whose deferred calls run
# first): this directory's tier labels, then the checks of every directory and of the list.
fepx_defer_test_tiers()
cmake_language(DEFER CALL fepx_check_quick_tests)

# Assertion module, linked into every test executable.
add_library(fepx_testing STATIC framework/assert_mod.f90)

# Mesh distribution in blocks on any number of ranks, for test programs (D108).
add_library(fepx_test_mesh STATIC framework/mesh_blocks_mod.f90)
target_link_libraries(fepx_test_mesh PUBLIC fepx_parallel)

# --- Assertion module self-tests -------------------------------------------------------------

# Checks that must pass.
fepx_add_unit_test(test_assert SOURCES framework/test_assert.f90)

# Checks that must fail: each case performs one passing and one failing check, and the output
# must report exactly one failure out of two checks.
add_executable(test_assert_failures framework/test_assert_failures.f90)
target_link_libraries(test_assert_failures PRIVATE fepx_testing)
set(assert_failure_cases
    true int32 int64 char
    scalar_exact scalar_abs scalar_rel scalar_nan scalar_inf_sign scalar_inf_finite
    rank1_value rank1_nan rank1_size rank2_value rank2_shape)
foreach(failure_case IN LISTS assert_failure_cases)
  add_test(NAME test_assert_failure_${failure_case}
           COMMAND test_assert_failures ${failure_case})
  set_tests_properties(test_assert_failure_${failure_case} PROPERTIES
                       PASS_REGULAR_EXPRESSION "test summary: 2 checks, 1 failed"
                       LABELS unit)
endforeach()

# A test that performs no check is reported as a failure.
add_test(NAME test_assert_failure_no_checks COMMAND test_assert_failures no_checks)
set_tests_properties(test_assert_failure_no_checks PROPERTIES
                     PASS_REGULAR_EXPRESSION "no checks were performed"
                     LABELS unit)

# A failing test exits with a non-zero status (PASS_REGULAR_EXPRESSION above ignores it).
add_test(NAME test_assert_failure_exit_status COMMAND test_assert_failures int32)
set_tests_properties(test_assert_failure_exit_status PROPERTIES
                     WILL_FAIL TRUE
                     LABELS unit)

# --- Core modules (src/core/) ----------------------------------------------------------------

fepx_add_unit_test(test_kinds SOURCES core/test_kinds.f90)
fepx_add_unit_test(test_constants SOURCES core/test_constants.f90)

# error_mod: each case of test_error_cases ends the program in a known way; check_run.cmake
# checks its exact exit status, its standard output and its standard error. Entries:
# name|status|message regex (both units)|extra standard-error regex after the message (or empty).
add_executable(test_error_cases core/test_error_cases.f90)
target_link_libraries(test_error_cases PRIVATE fepx_core)
set(error_cases
    "fatal_default|1|Error  : test_routine: expected failure|"
    "internal|1|Error  : test_routine: internal error [(]please report it[)]: a table out of order|"
    "fatal_handler|7|Error  : test_routine: failure through the handler|.*stopping handler called with status 1"
    "fatal_returning|1|Error  : test_routine: failure with a returning handler|.*returning handler called with status 1"
    "fatal_reset|1|Error  : test_routine: failure after reset|"
    "status_errmsg|1|Error  : test_routine: open file missing.cfg failed \\(status 5\\): No such file|"
    "status_no_errmsg|1|Error  : test_routine: read line 4 failed \\(status -3\\)[^:]|"
    "status_blank_errmsg|1|Error  : test_routine: read line 5 failed \\(status 2\\)[^:]|"
    "status_iostat|1|Error  : test_routine: open file this_file_does_not_exist.cfg failed \\(status [0-9-]+\\): [^ ]|")
foreach(error_case IN LISTS error_cases)
  string(REPLACE "|" ";" error_case "${error_case}")
  list(GET error_case 0 case_name)
  list(GET error_case 1 case_status)
  list(GET error_case 2 case_message)
  list(GET error_case 3 case_stderr_extra)
  add_test(NAME test_error_${case_name}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=$<TARGET_FILE:test_error_cases>;${case_name}"
                   -DEXPECTED_STATUS=${case_status}
                   "-DEXPECTED_STDOUT=${case_message}"
                   "-DEXPECTED_STDERR=${case_message}${case_stderr_extra}"
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_error_${case_name} PROPERTIES LABELS unit)
endforeach()

# check_status with a zero status returns and the program ends normally.
add_test(NAME test_error_status_ok
         COMMAND ${CMAKE_COMMAND}
                 "-DCOMMAND=$<TARGET_FILE:test_error_cases>;status_ok"
                 -DEXPECTED_STATUS=0
                 "-DEXPECTED_STDOUT=check_status returned"
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
set_tests_properties(test_error_status_ok PROPERTIES LABELS unit)

# log_mod: line format, rank selection, all-ranks tag.
fepx_add_unit_test(test_log SOURCES core/test_log.f90)

# --- Tensor algebra (src/math/) --------------------------------------------------------------

# Shared test helpers (random inputs, explicit 3x3/fourth-order algebra), compiled once: a source
# compiled into several executables would write the same module file concurrently.
add_library(fepx_test_math_utils STATIC math/tensor_test_utils_mod.f90)
target_link_libraries(fepx_test_math_utils PRIVATE fepx_core)

fepx_add_unit_test(test_tensor SOURCES math/test_tensor.f90)
target_link_libraries(test_tensor PRIVATE fepx_test_math_utils)
fepx_add_unit_test(test_rotation SOURCES math/test_rotation.f90)
target_link_libraries(test_rotation PRIVATE fepx_test_math_utils)
fepx_add_unit_test(test_linalg SOURCES math/test_linalg.f90)
target_link_libraries(test_linalg PRIVATE fepx_test_math_utils)
fepx_add_unit_test(test_time_integration SOURCES math/test_time_integration.f90)
target_link_libraries(test_time_integration PRIVATE fepx_test_math_utils)

# Orientation descriptors: reads the Neper fixtures (math/fixtures/neper_ori/README.md).
add_executable(test_orientation math/test_orientation.f90)
target_link_libraries(test_orientation PRIVATE fepx_testing fepx_core fepx_test_math_utils)
add_test(NAME test_orientation
         COMMAND test_orientation ${CMAKE_CURRENT_SOURCE_DIR}/math/fixtures/neper_ori)
set_tests_properties(test_orientation PROPERTIES LABELS unit)

# Orientation input validation: each case must end with a fatal error (status 1) and its message.
add_executable(test_orientation_errors math/test_orientation_errors.f90)
target_link_libraries(test_orientation_errors PRIVATE fepx_core)
set(orientation_error_cases
    "bad_descriptor|invalid orientation descriptor code 8"
    "bad_convention|invalid orientation convention code 3"
    "too_few_params|expected 3 parameters for euler-bunge, got 2"
    "bad_rotmat|rotmat orientation is not a rotation matrix"
    "reflection_rotmat|rotmat orientation is not a rotation matrix"
    "zero_axis|axis-angle orientation with a zero axis"
    "zero_quaternion|zero quaternion orientation"
    "bad_output|descriptor_from_rot_c2s: internal error [(]please report it[)]: invalid orientation descriptor code 0")
foreach(error_case IN LISTS orientation_error_cases)
  string(REPLACE "|" ";" error_case "${error_case}")
  list(GET error_case 0 case_name)
  list(GET error_case 1 case_message)
  add_test(NAME test_orientation_error_${case_name}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=$<TARGET_FILE:test_orientation_errors>;${case_name}"
                   -DEXPECTED_STATUS=1
                   "-DEXPECTED_STDOUT=Error  : .*${case_message}"
                   "-DEXPECTED_STDERR=Error  : .*${case_message}"
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_orientation_error_${case_name} PROPERTIES LABELS unit)
endforeach()

# --- Crystal module (src/material/) --------------------------------------------------------

# Slip systems and stiffness; compares with the reference Schmid tensors (material/fixtures/).
add_executable(test_crystal material/test_crystal.f90)
target_link_libraries(test_crystal PRIVATE fepx_testing fepx_core fepx_test_math_utils)
add_test(NAME test_crystal
         COMMAND test_crystal ${CMAKE_CURRENT_SOURCE_DIR}/material/fixtures/schmid_reference.txt)
set_tests_properties(test_crystal PROPERTIES LABELS unit)

# crystal_init input validation: fatal error (status 1) with its message.
add_executable(test_crystal_errors material/test_crystal_errors.f90)
target_link_libraries(test_crystal_errors PRIVATE fepx_core)
set(crystal_error_cases
    "hcp_no_c_over_a|internal error [(]please report it[)]: c_over_a must be given and positive for hcp"
    "bct_zero_c_over_a|internal error [(]please report it[)]: c_over_a must be given and positive for bct"
    "bad_structure|internal error [(]please report it[)]: invalid crystal structure code 5"
    "no_family|internal error [(]please report it[)]: no active slip family"
    "bad_family|internal error [(]please report it[)]: invalid slip family 2 for fcc"
    "family_twice|internal error [(]please report it[)]: slip family 1 given twice")
foreach(error_case IN LISTS crystal_error_cases)
  string(REPLACE "|" ";" error_case "${error_case}")
  list(GET error_case 0 case_name)
  list(GET error_case 1 case_message)
  add_test(NAME test_crystal_error_${case_name}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=$<TARGET_FILE:test_crystal_errors>;${case_name}"
                   -DEXPECTED_STATUS=1
                   "-DEXPECTED_STDOUT=Error  : crystal_init: ${case_message}"
                   "-DEXPECTED_STDERR=Error  : crystal_init: ${case_message}"
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_crystal_error_${case_name} PROPERTIES LABELS unit)
endforeach()

# Slip kinetics (step 6.1): cutoff, resolved shear stress, power law, [M*] and [M*_tan].
fepx_add_unit_test(test_slip_kinetics SOURCES material/test_slip_kinetics.f90)
target_link_libraries(test_slip_kinetics PRIVATE fepx_test_math_utils)

# Hardening laws and their time integration (step 6.4).
fepx_add_unit_test(test_hardening SOURCES material/test_hardening.f90)
target_link_libraries(test_hardening PRIVATE fepx_test_math_utils)
fepx_add_unit_test(test_hardening_km SOURCES material/test_hardening_km.f90)
target_link_libraries(test_hardening_km PRIVATE fepx_test_math_utils)
fepx_add_unit_test(test_hardening_irr SOURCES material/test_hardening_irr.f90)
target_link_libraries(test_hardening_irr PRIVATE fepx_test_math_utils)

# Crystal plasticity point solve (step 6.2).
fepx_add_unit_test(test_cp_model SOURCES material/test_cp_model.f90)
target_link_libraries(test_cp_model PRIVATE fepx_test_math_utils)
fepx_add_unit_test(test_cp_model_km SOURCES material/test_cp_model_km.f90)

fepx_add_unit_test(test_cp_model_irr SOURCES material/test_cp_model_irr.f90)
fepx_add_unit_test(test_cp_state_damping SOURCES material/test_cp_state_damping.f90)
# v2's yield-surface vertex tables and its increment-1 viscoplastic estimate (step 6.6).
fepx_add_unit_test(test_vertex_estimate SOURCES material/test_vertex_estimate.f90)
target_link_libraries(test_vertex_estimate PRIVATE fepx_test_math_utils)

# Computed yield-surface vertices (step 6.7); compares with KD25's implementation
# (material/fixtures/kd25_vertices.txt, README).
add_executable(test_yield_vertices material/test_yield_vertices.f90)
target_link_libraries(test_yield_vertices PRIVATE fepx_testing fepx_core fepx_test_math_utils)
add_test(NAME test_yield_vertices
         COMMAND test_yield_vertices ${CMAKE_CURRENT_SOURCE_DIR}/material/fixtures/kd25_vertices.txt)
set_tests_properties(test_yield_vertices PROPERTIES LABELS unit)

# Output quantities of the material state (step 6.5).
fepx_add_unit_test(test_material_outputs SOURCES material/test_material_outputs.f90)
target_link_libraries(test_material_outputs PRIVATE fepx_test_math_utils)

# Addendum_hardening §5 test 3 and addendum_point_verification (step 7.4a): the driver writes
# point-model histories for halved dt; the Python check integrates the rate equations
# independently (RK4) and compares, in four groups of cases that run in parallel: 1-6 (cube,
# no lattice spin), 7-10 (lattice spin, shear, saturation evolution), 11-14 (cyclic,
# precipitates, bcc, hcp), 15-18 (Kocks-Mecking limits).
add_executable(drive_hardening_point material/drive_hardening_point.f90)
target_link_libraries(drive_hardening_point PRIVATE fepx_core)
add_test(NAME hardening_integration_run
         COMMAND drive_hardening_point ${CMAKE_CURRENT_BINARY_DIR}/hardening_integration.txt)
set_tests_properties(hardening_integration_run PROPERTIES LABELS verification
                     FIXTURES_SETUP hardening_integration)
foreach(group IN ITEMS "check;1 2 3 4 5 6" "spin_check;7 8 9 10" "modifier_check;11 12 13 14"
                       "km_limit_check;15 16 17 18")
  list(GET group 0 name)
  list(GET group 1 cases)
  separate_arguments(cases)
  add_test(NAME hardening_integration_${name}
           COMMAND ${Python3_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/verification/check_hardening_integration.py
                   ${CMAKE_CURRENT_BINARY_DIR}/hardening_integration.txt --cases ${cases})
  set_tests_properties(hardening_integration_${name} PROPERTIES LABELS verification
                       FIXTURES_REQUIRED hardening_integration)
endforeach()

# Material-solver contract and the elastic-only model (step 4.1).
fepx_add_unit_test(test_elastic_model SOURCES material/test_elastic_model.f90)
target_link_libraries(test_elastic_model PRIVATE fepx_test_math_utils)

# --- Element library (src/element/) ---------------------------------------------------------

# 10-node tetrahedron and its quadrature rule; compares with the reference (element/fixtures/) and
# checks the Gmsh node map on the io fixtures and two Neper regression meshes.
add_executable(test_tet10 element/test_tet10.f90)
target_link_libraries(test_tet10 PRIVATE fepx_testing fepx_core fepx_test_math_utils)
add_test(NAME test_tet10
         COMMAND test_tet10 ${CMAKE_CURRENT_SOURCE_DIR}/element/fixtures
                 ${CMAKE_CURRENT_SOURCE_DIR}/io/fixtures
                 ${CMAKE_CURRENT_SOURCE_DIR}/regression/cases/meshes)
set_tests_properties(test_tet10 PROPERTIES LABELS unit)

# 20-node hexahedron and its 27-point rule; checks the Gmsh node order on a Neper hex mesh
# (element/fixtures/hex_cube.msh).
add_executable(test_hex20 element/test_hex20.f90)
target_link_libraries(test_hex20 PRIVATE fepx_testing fepx_core fepx_test_math_utils)
add_test(NAME test_hex20 COMMAND test_hex20 ${CMAKE_CURRENT_SOURCE_DIR}/element/fixtures)
set_tests_properties(test_hex20 PROPERTIES LABELS unit)

# Element descriptor and geometry kernels: reference fixture (element/fixtures/geometry_reference.txt),
# property and patch tests, and the Neper meshes.
add_executable(test_element element/test_element.f90)
target_link_libraries(test_element PRIVATE fepx_testing fepx_core fepx_test_math_utils)
add_test(NAME test_element
         COMMAND test_element ${CMAKE_CURRENT_SOURCE_DIR}/element/fixtures
                 ${CMAKE_CURRENT_SOURCE_DIR}/regression/cases/meshes)
set_tests_properties(test_element PROPERTIES LABELS unit)

# Element errors: each case must end with a fatal error (status 1) and its message.
add_executable(test_element_errors element/test_element_errors.f90)
target_link_libraries(test_element_errors PRIVATE fepx_core)
set(element_error_case_message
    "Error  : element_init: unsupported element type [(]Gmsh type code 4, supported: 11 and 17[)]")
add_test(NAME test_element_error_unknown_type
         COMMAND ${CMAKE_COMMAND}
                 "-DCOMMAND=$<TARGET_FILE:test_element_errors>;unknown_type"
                 -DEXPECTED_STATUS=1
                 "-DEXPECTED_STDOUT=${element_error_case_message}"
                 "-DEXPECTED_STDERR=${element_error_case_message}"
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
set_tests_properties(test_element_error_unknown_type PROPERTIES LABELS unit)

# --- Input files (src/io/) -------------------------------------------------------------------

fepx_add_unit_test(test_text_file SOURCES io/test_text_file.f90)

# Mesh parser on the fixtures (io/fixtures/README.md) and on two Neper regression meshes.
add_executable(test_mesh_file io/test_mesh_file.f90)
target_link_libraries(test_mesh_file PRIVATE fepx_testing fepx_core)
add_test(NAME test_mesh_file
         COMMAND test_mesh_file ${CMAKE_CURRENT_SOURCE_DIR}/io/fixtures
                 ${CMAKE_CURRENT_SOURCE_DIR}/regression/cases/meshes)
set_tests_properties(test_mesh_file PROPERTIES LABELS unit)

# Mesh parser errors: each case writes a modified copy of io/fixtures/tiny.msh in its own working
# directory and must end with a fatal error (status 1) and its message on both output units.
add_executable(test_mesh_file_errors io/test_mesh_file_errors.f90)
target_link_libraries(test_mesh_file_errors PRIVATE fepx_core)
set(mesh_error_cases
    "missing_file|Cannot open file .nonexistent.msh.: No such file"
    "missing_version|no [$]MeshVersion section: meshes written before Neper msh file version 2.4 are not supported"
    "version_2_2_3|mesh file version 2.2.3 is not supported: FEPX reads Neper msh file version 2.4 only"
    "version_2_3|mesh file version 2.3 is not supported"
    "format_binary|mesh format must be .2.2 0 8."
    "missing_format|no [$]MeshFormat section"
    "missing_elements|no [$]Elements section"
    "not_a_header|expected a section header"
    "unexpected_end|unexpected [$]EndNodes"
    "node_order|node ids must be 1, 2, [.][.][.] in order: expected 2, found 3"
    "node_values|node line: expected 4 values, found 3"
    "node_not_number|node line: expected numbers"
    "count_not_integer|number of nodes: expected integers"
    "node_count_footer|expected [$]EndNodes"
    "overflow_dofs|number of degrees of freedom [(]3 x nodes[)] 2400000000 exceeds the largest integer"
    "unused_node|node 15 belongs to no 3D element"
    "elt_first_order|element type 4 is not supported"
    "elt_tags|expected 3 tags"
    "elt_id_gap|3D element ids must be consecutive: expected 4, found 5"
    "elt_node_range|refers to a node id outside 1[.][.]14"
    "elt_none|no second-order 3D elements"
    "duplicate_section|section [$]Periodicity given twice"
    "partial_partition|some 3D elements have a partition and others not"
    "missing_node_parts|there is no [$]NodePartitions section"
    "node_part_range|partition 3 of node 5 outside 1[.][.]2"
    "nset_range|nset apex has a node id outside 1[.][.]14"
    "faset_wrong_elt|node 1 of face 1 is not a node of its element [(]id 4[)]"
    "faset_not_3d|refers to element id 2, which is not a 3D element"
    "periodicity_vector|periodicity vector components must be -1, 0 or 1"
    "ori_unknown_elset|elset 5 does not exist"
    "ori_twice|elset 3 given twice"
    "ori_count|has 1 lines for 2 elsets"
    "ori_descriptor|invalid orientation descriptor .euler:active."
    "ori_params|id 7 has 2 parameters, expected 3"
    "groups_entity|the group entity must be .elset."
    "truncated|unexpected end of file after line [0-9]+ [(]expected [$]EndGroups[)]"
    "groups_phase|groups [(]phases[)] must be positive"
    "aux_forbidden|section [$]Nodes is not allowed in an auxiliary file"
    "aux_crysym|[$]ElsetCrySym .cubic. differs from the mesh file's .triclinic."
    "aux_elt_ids|element 3 outside 1[.][.]2")
foreach(error_case IN LISTS mesh_error_cases)
  string(REPLACE "|" ";" error_case "${error_case}")
  list(GET error_case 0 case_name)
  list(GET error_case 1 case_message)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/io_errors/${case_name})
  file(MAKE_DIRECTORY ${case_dir})
  add_test(NAME test_mesh_file_error_${case_name}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=$<TARGET_FILE:test_mesh_file_errors>;${case_name};${CMAKE_CURRENT_SOURCE_DIR}/io/fixtures"
                   -DEXPECTED_STATUS=1
                   "-DEXPECTED_STDOUT=Error  : .*${case_message}"
                   "-DEXPECTED_STDERR=Error  : .*${case_message}"
                   -DWORKING_DIR=${case_dir}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_mesh_file_error_${case_name} PROPERTIES LABELS unit)
endforeach()

# Configuration reader: keyword table consistency and typed configurations of valid files.
fepx_add_unit_test(test_config SOURCES io/test_config.f90)

# The option-compatibility rules (cfg_rules_mod, step 9.9, D104).
fepx_add_unit_test(test_cfg_rules SOURCES io/test_cfg_rules.f90)

# The manual's generated keyword pages (cfg_manual_mod, step 9.9, D104).
fepx_add_unit_test(test_cfg_manual SOURCES io/test_cfg_manual.f90)

# Configuration reader errors (dev/cfg_keywords.md §1a): each case writes the program's valid
# base file with one change to case.cfg in its own working directory and must end with a fatal
# error (status 1) and its message on both output units. Entries: name|message regex.
add_executable(test_config_errors io/test_config_errors.f90)
target_link_libraries(test_config_errors PRIVATE fepx_core)
set(config_error_cases
    "unknown_keyword|case.cfg:20: unknown keyword .foo."
    "planned_keyword|case.cfg:20: temperature is not available in this version of FEPX"
    "output_format_invalid|case.cfg:20: output_format: invalid value .vtk. [(]expected one of: ascii hdf5[)]"
    "output_format_twice|case.cfg:20: output_format: give ascii, hdf5, or both once each"
    "output_format_three|case.cfg:20: output_format: give ascii, hdf5, or both once each"
    "load_number_of_incrs|case.cfg:17: step 2: number_of_incrs cannot be used with a load or stress target at a strain rate, or a strain target at a load or stress rate [(]the duration of the step is not known: give dtime or dstrain [(]the largest increment[)] and dtime_min. this step: a load target[)]"
    "load_no_dtime_min|case.cfg: step 2: a load target needs dtime_min, the smallest time increment of the step"
    "load_dtime_min_dash|case.cfg: step 2: a stress target needs dtime_min, the smallest time increment of the step [(]dtime_min, line 22, has .-.[)]"
    "dtime_min_not_load|case.cfg:20: step 1: dtime_min is used only with a load or stress target at a strain rate, or a strain target at a load or stress rate [(]give .-. for the other steps[)]"
    "dtime_min_large|case.cfg:22: step 2: dtime_min 0.500000 s is larger than dtime 0.100000 s"
    "two_targets_load|case.cfg: step 1, component 33: two targets [(]target_strain 33, line 16, and target_load 33, line 20[)]"
    "load_other_both|case.cfg: step 1, component 11: strain_rate 11 [(]line 21[)] and target_load 11 [(]line 22[)] both given. a component other than the control component 33 takes one, .-. for the others"
    "load_other_zero_rate|case.cfg: step 1, component 11: target_load 11 with strain_rate 33 0: the load-controlled face velocities are found by perturbations of min_pert_frac times the control face velocity, which must not be 0"
    "min_pert_frac_zero|case.cfg:22: min_pert_frac: value 0 out of range [(]expected > 0[)]"
    "load_other_target_strain|case.cfg: step 2, component 11: target_strain 11 [(]line 21[)] in a step that ends at target_load 33: the duration of the step is not known. give strain_rate 11"
    "dtime_factor_no_load|case.cfg:20: dtime_factor is used only with a target_load or target_stress"
    "dtime_factor_low|case.cfg:20: dtime_factor: value 0.9 out of range [(]expected >= 1[)]"
    "load_tol_no_load|case.cfg:20: load_tol is used only with a target_load or load_rate row"
    "stress_tol_no_stress|case.cfg:20: stress_tol is used only with a target_stress or stress_rate row"
    "load_tol_zero|case.cfg:23: load_tol: value 0 out of range [(]expected > 0[)]"
    "stress_tol_negative|case.cfg:23: stress_tol: value -0.1 out of range [(]expected > 0[)]"
    "max_bc_iter_no_load|case.cfg:20: max_bc_iter is used only with a target_load or target_stress"
    "max_bc_iter_zero|case.cfg:23: max_bc_iter: value 0 out of range [(]expected >= 1[)]"
    "rate_two|case.cfg: step 1, component 33: two rates [(]strain_rate 33, line 15, and load_rate 33, line 20[)]. give one, .-. for the other"
    "rate_general|case.cfg:20: load_rate cannot be used with boundary_conditions general"
    "load_rate_target_stress|case.cfg: step 1, component 33: load_rate 33 [(]line 20[)] cannot be used with a target_stress of the same component in the same step [(]a load_rate ends at a target_load, target_strain, target_time or target_duration. target_stress 33, line 21[)]"
    "stress_rate_target_load|case.cfg: step 1, component 33: stress_rate 33 [(]line 20[)] cannot be used with a target_load of the same component in the same step [(]a stress_rate ends at a target_stress, target_strain, target_time or target_duration. target_load 33, line 21[)]"
    "dwell_target_load|case.cfg: step 2, component 33: load_rate 33 0 [(]line 20[)] never reaches target_load 33: a rate 0 holds the value [(]a dwell[)]. give target_time or target_duration"
    "dwell_first|case.cfg: step 1, component 33: a dwell [(]load_rate 33 0, line 20[)] cannot be the first step"
    "dwell_number_of_incrs|case.cfg:17: step 2: number_of_incrs cannot be used with a dwell [(]a step with load_rate 0 or stress_rate 0 on the control component: give dtime, its largest increment. its increments follow max_strain_incr[)]"
    "rate_dstrain|case.cfg:22: step 1: dstrain is used only with a strain_rate on the control component [(]not in a step at a load_rate or stress_rate: give number_of_incrs or dtime[)]"
    "open_strain_no_dtime_min|case.cfg: step 1: a strain target at a load rate needs dtime_min"
    "open_strain_other_load|case.cfg: step 1, component 11: target_load 11 [(]line 23[)] in a step that ends at target_strain 33 at a load_rate 33: the duration of the step is not known"
    "duration_two_targets|case.cfg: step 1, component 33: two targets [(]target_strain 33, line 16, and target_duration, line 20[)]"
    "duration_zero|case.cfg:20: target_duration: value 0.0 out of range [(]expected > 0[)]"
    "max_strain_incr_no_dwell|case.cfg:20: max_strain_incr is used only with a dwell"
    "max_strain_incr_zero|case.cfg:20: max_strain_incr: value 0 out of range [(]expected > 0[)]"
    "min_pert_frac_no_load_rows|case.cfg:20: min_pert_frac is used only with a target_load or target_stress row on a component other than the control component, or a load_rate or stress_rate row"
    "check_necking_invalid|case.cfg:20: check_necking: invalid value .yes. [(]expected one of: on off[)]"
    "check_necking_no_component|case.cfg:20: check_necking on is used only with a loaded component [(]a target_strain, target_load or target_stress row. the load of its face x1, y1 or z1 is checked[)]"
    "planned_value|case.cfg:20: print strain_th is not available in this version of FEPX"
    "slip_with_elastic|case.cfg:20: m is used only with model crystal_plasticity [(]phase 1 is elastic[)]"
    "cp_no_m|case.cfg: phase 1 [(]hcp, line 3[)]: m is required with model crystal_plasticity"
    "cp_no_gammadot_0|case.cfg: phase 1 [(]hcp, line 3[)]: gammadot_0 is required with model crystal_plasticity"
    "cp_no_c_over_a|case.cfg: phase 1 [(]hcp, line 3[)]: c_over_a is required with model crystal_plasticity"
    "cp_m_count|case.cfg:20: m: 2 value[(]s[)]. give one per active slip family [(]3[)] or one per slip system [(]18[)]"
    "cp_m_single_several_families|case.cfg:20: m: 1 value[(]s[)]. give one per active slip family [(]3[)] or one per slip system [(]18[)]"
    "cp_gammadot_0_count_subset|case.cfg:22: gammadot_0: 2 value[(]s[)]. give one per active slip family [(]1[)] or one per slip system [(]3[)]"
    "cp_family_unknown|case.cfg:20: slip_families: .[{]111[}]<110>. is not a slip family of crystal_type hcp [(]expected: [{]0001[}]<11-20>, [{]10-10[}]<11-20>, [{]10-11[}]<11-23>, or all[)]"
    "cp_family_twice|case.cfg:20: slip_families: [{]0001[}]<11-20> given twice"
    "cp_all_combined|case.cfg:20: slip_families: all cannot be combined with other families"
    "cp_family_case|case.cfg:20: slip_families: .ALL. is not a slip family of crystal_type hcp"
    "cp_m_zero|case.cfg:20: m: value 0.0 out of range [(]expected > 0[)]"
    "cp_gammadot_0_negative|case.cfg:20: gammadot_0: value -1.0 out of range [(]expected > 0[)]"
    "list_empty|case.cfg:20: slip_families: value[(]s[)] missing"
    "cp_no_tau_0|case.cfg: phase 1 [(]hcp, line 3[)]: tau_0 is required with model crystal_plasticity"
    "cp_tau_0_count|case.cfg:22: tau_0: 2 value[(]s[)]. give one per active slip family [(]3[)] or one per slip system [(]18[)]"
    "cp_tau_0_zero|case.cfg:20: tau_0: value 0 out of range [(]expected > 0[)]"
    "tau_0_with_elastic|case.cfg:20: tau_0 is used only with model crystal_plasticity [(]phase 1 is elastic[)]"
    "initialization_invalid|case.cfg:20: initialization: invalid value .plastic. [(]expected one of: elastic_step isotropic_viscoplastic[)]"
    "sx_tol_zero|case.cfg:20: point_tol: value 0.0 out of range [(]expected > 0[)]"
    "tol_state_zero|case.cfg:20: point_tol_state: value 0 out of range [(]expected > 0[)]"
    "tol_hardening_negative|case.cfg:20: point_tol_hardening: value -1e-4 out of range [(]expected > 0[)]"
    "sx_iters_zero|case.cfg:20: point_max_iters_newton: value 0 out of range [(]expected >= 1[)]"
    "hard_with_elastic|case.cfg:20: h_0 is used only with model crystal_plasticity [(]phase 1 is elastic[)]"
    "hard_missing_h_0|case.cfg: phase 1 [(]hcp, line 3[)]: h_0 is required with hardening saturation_isotropic"
    "hard_two_values|case.cfg:23: h_0: hardening saturation_isotropic takes one value [(]one hardening increment for all systems, set by slip system 1[)]"
    "hard_tau_s_low|case.cfg:24: tau_s must be greater than tau_0 of slip system 1 [(]390.000[)]"
    "hard_none_with_h_0|case.cfg:24: h_0 is used only with a saturation hardening law [(]hardening saturation_isotropic or saturation_anisotropic. phase 1 has hardening none[)]"
    "mod_with_none|case.cfg:24: saturation_evolution is used only with a saturation hardening law [(]hardening saturation_isotropic or saturation_anisotropic. phase 1[)]"
    "mod_with_elastic|case.cfg:20: cyclic_hardening is used only with a saturation hardening law [(]hardening saturation_isotropic or saturation_anisotropic. phase 1[)]"
    "sat_evo_missing|case.cfg: phase 1 [(]hcp, line 3[)]: m_prime is required with saturation_evolution on"
    "sat_evo_off_param|case.cfg:27: m_prime is used only with saturation_evolution on"
    "cyclic_off_param|case.cfg:27: cyclic_a is used only with cyclic_hardening on"
    "precipitate_off_param|case.cfg:27: a_p is used only with precipitate_strengthening on"
    "sat_evo_tau_s_values|case.cfg:25: tau_s: saturation_evolution on takes one value"
    "cyclic_a_range|case.cfg:28: cyclic_a must lie between 0 and 1"
    "cyclic_missing|case.cfg: phase 1 [(]hcp, line 3[)]: cyclic_c is required with cyclic_hardening on"
    "precip_incomplete|case.cfg: phase 1 [(]hcp, line 3[)]: precipitate_strengthening on needs the cutting parameters"
    "precip_spacing|case.cfg:29: precipitate bowing needs a spacing r_p sqrt[(]pi/f_p[)] larger than 2 r_p"
    "precip_tau_s_low|case.cfg: phase 1 [(]hcp, line 3[)]: tau_s must be greater than tau_0 [+] the precipitate strength"
    "km_missing|case.cfg: phase 1 [(]hcp, line 3[)]: km_n is required with hardening kocks_mecking"
    "km_with_saturation|case.cfg:27: km_b is used only with hardening kocks_mecking [(]phase 1[)]"
    "km_with_h_0|case.cfg:31: h_0 is used only with a saturation hardening law [(]hardening saturation_isotropic or saturation_anisotropic. phase 1 has hardening kocks_mecking[)]"
    "km_with_matrix|case.cfg:31: interaction_matrix is used only with a saturation hardening law [(]hardening saturation_isotropic or saturation_anisotropic. phase 1 has hardening kocks_mecking[)]"
    "km_modifier_off|case.cfg:31: precipitate_strengthening is used only with a saturation hardening law [(]hardening saturation_isotropic or saturation_anisotropic. phase 1 has hardening kocks_mecking, which takes no modifiers[)]"
    "km_modifier_on|case.cfg:31: cyclic_hardening is used only with a saturation hardening law [(]hardening saturation_isotropic or saturation_anisotropic. phase 1 has hardening kocks_mecking, which takes no modifiers[)]"
    "km_rho_0_zero|case.cfg:26: km_rho_0: value 0 out of range [(]expected > 0[)]"
    "km_n_zero|case.cfg:30: km_n: value 0 out of range [(]expected > 0[)]"
    "km_h_n_negative|case.cfg:25: km_h_n: value -0.1 out of range [(]expected >= 0[)]"
    "km_print_without_km|case.cfg:27: print km_rho is used only with a phase with hardening kocks_mecking"
    "irr_print_without_irr|case.cfg:31: print dmgrate is used only with a phase with irradiation on"
    "irr_missing|case.cfg: phase 1 [(]bcc, line 3[)]: irr_loop_normals is required with irradiation on"
    "irr_without_km|case.cfg:27: irradiation is used only with hardening kocks_mecking [(]phase 1[)]"
    "irr_param_off|case.cfg:31: irr_eta is used only with irradiation on [(]phase 1[)]"
    "irr_hcp|case.cfg:31: irradiation on is used only with crystal_type fcc and bcc [(]the dislocation loop normals are defined only for fcc and bcc. phase 1 has crystal_type hcp[)]"
    "irr_normals_bad|case.cfg:33: irr_loop_normals: .<112>. is not a loop-normal family [(]expected: <111>, <100>, <110>[)]"
    "irr_normals_twice|case.cfg:33: irr_loop_normals: <111> given twice"
    "irr_diameter_zero|case.cfg:32: irr_loop_diameter: value 0 out of range [(]expected > 0[)]"
    "irr_density_negative|case.cfg:31: irr_loop_density: value -1 out of range [(]expected >= 0[)]"
    "hard_iso_per_family|case.cfg:24: tau_s: hardening saturation_isotropic takes one value"
    "hard_iso_matrix|case.cfg:26: interaction_matrix is used only with hardening saturation_anisotropic"
    "hard_aniso_tau_s_low|case.cfg:25: tau_s must be greater than tau_0 for every slip system [(]system 7: tau_s 500.000, tau_0 663.000[)]"
    "hard_aniso_count|case.cfg:24: h_0: 2 value[(]s[)]. give one value, one per active slip family [(]3[)] or one per slip system [(]18[)]"
    "matrix_rows|case.cfg:27: interaction_matrix: 1 rows. give 18 [(]one per slip system[)]"
    "matrix_values|case.cfg:27: interaction_matrix: row 1 has 3 values. give 18"
    "matrix_negative|case.cfg:27: interaction_matrix: entry [(]1, 1[)] is negative. entries must be >= 0"
    "matrix_not_number|case.cfg:27: interaction_matrix: .x. is not a number"
    "subincr_depth_off|case.cfg:20: subincrement_max_depth is used only with subincrementation on"
    "subincr_successive_off|case.cfg:21: subincrement_max_successive is used only with subincrementation on"
    "subincr_depth_zero|case.cfg:21: subincrement_max_depth: value 0 out of range [(]expected >= 1[)]"
    "subincr_invalid|case.cfg:20: subincrementation: invalid value .yes. [(]expected one of: on off[)]"
    "nr_tol_switch_negative|case.cfg:20: global_nr_switch: value -1e-2 out of range [(]expected >= 0[)]"
    "nr_tol_conv_zero|case.cfg:20: global_nr_slow_ratio: value 0 out of range [(]expected > 0[)]"
    "planned_print|case.cfg:20: print temp is not available in this version of FEPX"
    "invalid_word|case.cfg:4: crystal_type: invalid value .xyz. [(]expected one of: fcc bcc hcp bct[)]"
    "phase_scope|case.cfg:1: c11 must follow a .phase. line"
    "duplicate|case.cfg:21: cg_tol given twice [(]lines 20 and 21[)]"
    "duplicate_row|case.cfg:20: strain_rate 33 given twice [(]lines 15 and 20[)]"
    "duplicate_phase_keyword|case.cfg:7: c11 given twice [(]lines 6 and 7[)]"
    "phase_twice|case.cfg:20: phase 1 defined twice [(]lines 3 and 20[)]"
    "value_count|case.cfg:20: cg_tol: expected 1 value[(]s[)], found 2"
    "value_missing|case.cfg:20: cg_tol: expected 1 value[(]s[)], found 0"
    "flag_value|case.cfg:20: read_ori_from_file takes no value"
    "repeat_count|case.cfg:20: print: expected 1 value[(]s[)]"
    "not_integer|case.cfg:20: cg_max_iters: .1.5. is not an integer"
    "not_number|case.cfg:6: c11: .abc. is not a number"
    "out_of_range|case.cfg:6: c11: value -5.0 out of range [(]expected > 0[)]"
    "out_of_range_open|case.cfg:20: cg_tol: value 0.0 out of range [(]expected > 0 and <= 1[)]"
    "dash_not_allowed|case.cfg:20: print_results: .-. is not allowed"
    "row_length|case.cfg:15: strain_rate: 1 values for 2 steps"
    "no_steps|case.cfg:14: strain_rate: per-step row without number_of_steps"
    "component_invalid|case.cfg:15: strain_rate: invalid component .44. [(]expected 11, 22, 33, or with periodic conditions 23, 13, 12, 32, 31, 21[)]"
    "component_shear|case.cfg:15: strain_rate 12: shear components are used only with boundary_conditions periodic"
    "control_shear|case.cfg:20: control_component 12: shear components are used only with boundary_conditions periodic"
    "component_missing|case.cfg:15: strain_rate: expected a component"
    "long_token|case.cfg:18: value .a+. is longer than 64 characters"
    "missing_required|case.cfg: boundary_conditions is required"
    "missing_phase_keyword|case.cfg: phase 1 [(]line 3[)]: model is required"
    "phase_range|case.cfg:20: phase 2 is outside 1[.][.]1 [(]number_of_phases[)]"
    "phase_missing|case.cfg: phase 2 is not defined [(]number_of_phases 2[)]"
    "missing_constant|case.cfg: phase 1 [(]hcp, line 3[)]: c33 is required"
    "extra_constant|case.cfg:10: c66 is not a constant of crystal_type hcp"
    "c_over_a_cubic|case.cfg:9: c_over_a is used only with crystal_type hcp and bct"
    "not_positive_definite|case.cfg: phase 1 [(]line 3[)]: the elastic constants are not positive definite"
    "two_components|case.cfg:20: target_strain 11: uniaxial_minimal loads one component, and line 15 gives 33"
    "no_strain_rate|case.cfg: step 1, component 33: no rate [(]boundary_conditions uniaxial_minimal needs a strain_rate 33 row[)]"
    "step_no_rate|case.cfg: step 2, component 33: no rate [(]strain_rate 33, line 15, has .-.[)]"
    "step_no_target|case.cfg: step 2, component 33: no target [(]target_strain 33, target_load 33, target_stress 33, target_time or target_duration[)]"
    "step_two_targets|case.cfg: step 2, component 33: two targets [(]target_strain 33, line 16, and target_time, line 20[)]"
    "general_no_target|case.cfg: step 1: no target [(]target_strain <component>, target_load <component>, target_stress <component>, target_time or target_duration[)]"
    "general_two_components|case.cfg:20: target_strain 11: general loads one component, and line 15 gives 33"
    "triaxial_no_control|case.cfg: boundary_conditions triaxial needs control_component <11, 22 or 33>"
    "control_no_rows|case.cfg:20: control_component 11 has no strain_rate or target_strain rows"
    "control_mismatch|case.cfg:20: control_component 11 is not the loaded component 33 of uniaxial_minimal"
    "control_general|case.cfg:20: control_component cannot be used with boundary_conditions general"
    "triaxial_both|case.cfg: step 1, component 11: strain_rate 11 [(]line 21[)] and target_strain 11 [(]line 22[)] both given"
    "triaxial_neither|case.cfg: step 2, component 11: no strain_rate, target_strain, target_load, target_stress, load_rate or stress_rate"
    "incr_none|case.cfg: step 2: no increments [(]number_of_incrs, dtime or dstrain[)]"
    "incr_two|case.cfg:20: step 1: number_of_incrs [(]line 17[)] and dtime both give the increments"
    "print_results_value|case.cfg:20: print_results: value 2 out of range"
    "set_bc_with_preset|case.cfg:20: set_bc cannot be used with boundary_conditions uniaxial_grip, uniaxial_symmetry, uniaxial_minimal or triaxial [(]the preset sets the velocities. given: boundary_conditions uniaxial_minimal[)]"
    "general_rate|case.cfg:15: strain_rate cannot be used with boundary_conditions general"
    "general_no_set_bc|case.cfg: boundary_conditions general needs set_bc vel lines"
    "periodic_kinematics_not_periodic|case.cfg:20: periodic_kinematics is used only with boundary_conditions periodic"
    "periodic_precond_not_periodic|case.cfg:20: periodic_preconditioner is used only with boundary_conditions periodic"
    "periodic_no_control|case.cfg: boundary_conditions periodic needs control_component <ij>"
    "periodic_control_not_imposed|case.cfg:20: control_component 22 has no strain_rate row"
    "periodic_small_transposed|case.cfg:21: strain_rate 21 cannot be used with periodic_kinematics small_strain [(]it takes the shear components 23, 13 and 12, one per pair. deformation_gradient_rate takes the nine components[)]"
    "periodic_rate_dash|case.cfg:21: strain_rate 12: step 2 has .-.: with periodic conditions a component is imposed"
    "periodic_target_other|case.cfg:22: target_strain 11: with boundary_conditions periodic only the control component 33 has strain targets"
    "periodic_load|case.cfg:21: target_load cannot be used with boundary_conditions periodic"
    "periodic_mpc|case.cfg:21: multi-point constraints cannot be used with boundary_conditions periodic"
    "cg_precond_value|case.cfg:20: cg_preconditioner: invalid value .ilu. [(]expected one of: jacobi block_jacobi none[)]"
    "cg_precond_periodic|case.cfg:21: cg_preconditioner block_jacobi cannot be used with boundary_conditions periodic [(]only jacobi[)]"
    "cg_precond_mpc|case.cfg:21: cg_preconditioner none cannot be used with multi-point constraints [(]set_mpc, set_mpc1: only jacobi[)]"
    "periodic_necking|case.cfg:21: check_necking on cannot be used with boundary_conditions periodic [(]no face loads[)]"
    "periodic_two_targets|case.cfg: step 1, component 33: two targets [(]target_strain 33, line 16, and target_time, line 21[)]"
    "set_bc_kind|case.cfg:20: set_bc: invalid kind .speed. [(]expected vel[)]"
    "set_bc_pairs|case.cfg:20: set_bc vel: expected <nset> and one to three"
    "set_bc_direction|case.cfg:20: set_bc vel: invalid direction .w. [(]expected x, y or z[)]"
    "set_bc_value|case.cfg:20: set_bc vel: .fast. is not a number"
    "set_bc_twice|case.cfg:20: set_bc vel: x0 x given twice [(]lines 19 and 20[)]"
    "mpc_no_vel|case.cfg:20: set_mpc needs vel: set_mpc vel <nset> <x.y.z> <coefficient> <offset>"
    "mpc1_kind|case.cfg:20: set_mpc1: invalid kind .speed. [(]expected vel[)]"
    "mpc1_count|case.cfg:20: set_mpc1 vel: expected <nset> and one to three directions"
    "mpc_count|case.cfg:20: set_mpc vel: expected <nset> and one to three <x.y.z> <coefficient> <offset> triples"
    "mpc_direction|case.cfg:20: set_mpc1 vel: invalid direction .w. [(]expected x, y or z[)]"
    "mpc_zero_coeff|case.cfg:20: set_mpc vel: the coefficient must be nonzero"
    "mpc_value|case.cfg:20: set_mpc vel: .fast. is not a number"
    "mpc_twice|case.cfg:21: set_mpc vel: x1 x given twice [(]lines 20 and 21[)]"
    "mpc_triaxial|case.cfg:21: multi-point constraints cannot be used with boundary_conditions triaxial"
    "mpc_load_target|case.cfg:21: multi-point constraints cannot be used with load or stress targets or rates [(]target_load, target_stress, load_rate, stress_rate: load control[)]"
    "print_twice|case.cfg:20: print coo given twice"
    "restart_invalid|case.cfg:20: restart: invalid value .yes. [(]expected one of: on off[)]")
# Output format hdf5 in a build without the HDF5 output (step 9.1).
if(NOT FEPX_HAVE_HDF5)
  list(APPEND config_error_cases
       "output_format_no_hdf5|case.cfg:20: output_format hdf5: this FEPX build has no HDF5 output [(]configure with -DFEPX_ENABLE_HDF5=ON[)]")
endif()
foreach(error_case IN LISTS config_error_cases)
  string(REPLACE "|" ";" error_case "${error_case}")
  list(GET error_case 0 case_name)
  list(GET error_case 1 case_message)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/config_errors/${case_name})
  file(MAKE_DIRECTORY ${case_dir})
  add_test(NAME test_config_error_${case_name}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=$<TARGET_FILE:test_config_errors>;${case_name}"
                   -DEXPECTED_STATUS=1
                   "-DEXPECTED_STDOUT=Error  : .*${case_message}"
                   "-DEXPECTED_STDERR=Error  : .*${case_message}"
                   -DWORKING_DIR=${case_dir}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_config_error_${case_name} PROPERTIES LABELS unit)
endforeach()
# The unmodified base file must be accepted.
set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/config_errors/valid)
file(MAKE_DIRECTORY ${case_dir})
add_test(NAME test_config_error_valid
         COMMAND ${CMAKE_COMMAND}
                 "-DCOMMAND=$<TARGET_FILE:test_config_errors>;valid"
                 -DEXPECTED_STATUS=0
                 -DWORKING_DIR=${case_dir}
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
set_tests_properties(test_config_error_valid PROPERTIES LABELS unit)
# A file that only warns: print forces with periodic conditions (step 8.3b, Q-P11; the rule
# print_forces_periodic, step 9.9).
set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/config_errors/periodic_print_forces)
file(MAKE_DIRECTORY ${case_dir})
add_test(NAME test_config_warning_periodic_print_forces
         COMMAND ${CMAKE_COMMAND}
                 "-DCOMMAND=$<TARGET_FILE:test_config_errors>;periodic_print_forces"
                 -DEXPECTED_STATUS=0
                 "-DEXPECTED_STDOUT=Warning: - case.cfg:21: print forces is ignored with boundary_conditions periodic [(]no face forces file. the macroscopic stress is written to results/macro[)]"
                 -DWORKING_DIR=${case_dir}
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
set_tests_properties(test_config_warning_periodic_print_forces PROPERTIES LABELS unit)

# --- MPI -------------------------------------------------------------------------------------

# Environment, reductions, broadcast and rank-0 logging through the wrappers.
fepx_add_mpi_test(test_parallel SOURCES parallel/test_parallel.f90 RANKS 1 2 4)

# Mesh distribution (also par_scatterv): partition table, not sequential too (D108), and blocks;
# 3 ranks also gives ranks without elements on tiny.msh.
fepx_add_mpi_test(test_mesh_distribute SOURCES parallel/test_mesh_distribute.f90 RANKS 1 2 3 4
                   ARGS ${CMAKE_CURRENT_SOURCE_DIR}/io/fixtures
                        ${CMAKE_CURRENT_SOURCE_DIR}/regression/cases/meshes)

# Local numbering, maps and communication graph (also par_alltoall/par_alltoallv).
fepx_add_mpi_test(test_mesh_local SOURCES parallel/test_mesh_local.f90 RANKS 1 2 3 4
                   ARGS ${CMAKE_CURRENT_SOURCE_DIR}/io/fixtures
                        ${CMAKE_CURRENT_SOURCE_DIR}/regression/cases/meshes)

# Halo exchange (forward, reverse-add) and distributed dot/norm: exact agreement with the serial
# assembly on every rank count.
fepx_add_mpi_test(test_halo SOURCES parallel/test_halo.f90 RANKS 1 2 3 4
                   ARGS ${CMAKE_CURRENT_SOURCE_DIR}/io/fixtures
                        ${CMAKE_CURRENT_SOURCE_DIR}/regression/cases/meshes)

# Error and abort cases: each run must end with an error, on every rank, without hanging (CTest
# timeout), with the message on standard output and standard error (check_run.cmake). Status
# "abort": the run must end through MPI_Abort (the abort handler of par_init) with a non-zero
# status (its value is up to the MPI implementation). (Until step 9.7 the handler wrote a line
# "MPI_Abort: ending all ranks", which these tests checked; it was removed with the terminal
# output of D94, and the launchers' own reports differ, so only the status and the absence of a
# hang are checked.) A message from a rank other than 0 is tagged "[Rank k]" (error_mod, Q62):
# the cases where only the last rank stops (one_rank, collective) on more than one rank require
# the tag on standard output, which shows that the message of a rank other than 0 reaches the
# log. "nonzero": any failing status, no MPI_Abort (after par_finalize). Every case that starts MPI runs through the launcher:
# singleton startup (an MPI program run without the launcher) hangs with some installations
# (seen with Open MPI 5.0.10). not_started does not start MPI and
# runs directly, with error_mod's exact status 1.
# Entries: name|rank counts (empty: run without the launcher)|expected status|message.
add_executable(test_parallel_errors parallel/test_parallel_errors.f90)
target_link_libraries(test_parallel_errors PRIVATE fepx_parallel)
set(parallel_error_cases
    "one_rank|1 2 4|abort|test_routine: error on one rank"
    "all_ranks|1 2 4|abort|test_routine: error on every rank"
    "collective|2 4|abort|test_routine: error during a reduction"
    "size_mismatch|1 2|abort|par_sum: internal error [(]please report it[)]: array sizes differ: contribution size 3, result size 2"
    "double_init|1|abort|par_init: internal error [(]please report it[)]: MPI is already initialized"
    "not_started||1|par_sum: internal error [(]please report it[)]: MPI is not started \\(call par_init first\\)"
    "after_finalize|1|nonzero|test_routine: error after finalize"
    "halo_size|1 2|abort|halo_forward: internal error [(]please report it[)]: vector size 5, expected 3 x 2"
    "dist_dot_size|1|abort|dist_dot: internal error [(]please report it[)]: num_owned exceeds the vector size")
foreach(error_case IN LISTS parallel_error_cases)
  string(REPLACE "|" ";" error_case "${error_case}")
  list(GET error_case 0 case_name)
  list(GET error_case 1 case_ranks)
  list(GET error_case 2 case_status)
  list(GET error_case 3 case_message)
  if(case_ranks STREQUAL "")
    set(case_runs "serial")
  else()
    string(REPLACE " " ";" case_runs "${case_ranks}")
  endif()
  foreach(case_run IN LISTS case_runs)
    if(case_run STREQUAL "serial")
      set(test_name test_parallel_error_${case_name})
      set(command $<TARGET_FILE:test_parallel_errors> ${case_name})
      set(test_processors 1)
    else()
      set(test_name test_parallel_error_${case_name}_np${case_run})
      fepx_mpiexec_command(command ${case_run} $<TARGET_FILE:test_parallel_errors>)
      list(APPEND command ${case_name})
      set(test_processors ${case_run})
    endif()
    set(expected_message "Error  : (\\[Rank [0-9]+\\] )?${case_message}")
    if((case_name STREQUAL "one_rank" OR case_name STREQUAL "collective") AND
       NOT case_run STREQUAL "serial" AND case_run GREATER 1)
      set(expected_message "Error  : \\[Rank [1-9][0-9]*\\] ${case_message}")
    endif()
    set(expected_status ${case_status})
    if(case_status STREQUAL "abort")
      set(expected_status nonzero)
    endif()
    add_test(NAME ${test_name}
             COMMAND ${CMAKE_COMMAND}
                     "-DCOMMAND=${command}"
                     -DEXPECTED_STATUS=${expected_status}
                     "-DEXPECTED_STDOUT=${expected_message}"
                     "-DEXPECTED_STDERR=${expected_message}"
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
    set_tests_properties(${test_name} PROPERTIES
                         TIMEOUT ${fepx_timeout_60}
                         PROCESSORS ${test_processors}
                         LABELS mpi)
  endforeach()
endforeach()

# --- Finite element operations (src/fe/) -----------------------------------------------------

# Element kinematics, element matrix and vector (serial).
fepx_add_unit_test(test_element_matrix SOURCES fe/test_element_matrix.f90)
target_link_libraries(test_element_matrix PRIVATE fepx_test_math_utils)

# EBE matrix-vector product, Jacobi-PCG and Dirichlet presets on 1, 2 and 4 ranks.
fepx_add_mpi_test(test_fe_solve SOURCES fe/test_fe_solve.f90 RANKS 1 2 4
                   ARGS ${CMAKE_CURRENT_SOURCE_DIR}/regression/cases/meshes)

# CG preconditioners (step 9.2): the nodal blocks against a serial assembly (tet10, hex20), the
# preconditioner step, invariance under nodal rotations, solves, a block not positive definite.
fepx_add_mpi_test(test_pcg_precond SOURCES fe/test_pcg_precond.f90 RANKS 1 2 4
                   ARGS ${CMAKE_CURRENT_SOURCE_DIR}/regression/cases/meshes
                        ${CMAKE_CURRENT_SOURCE_DIR}/verification/cases/meshes)

# Velocity constraints (step 8.3a): the MPC table, T and T^T, the reduced operator and its exact
# diagonal against a serial assembly, the propagation of velocity conditions, constrained solves.
fepx_add_mpi_test(test_constraint SOURCES fe/test_constraint.f90 RANKS 1 2 4
                   ARGS ${CMAKE_CURRENT_SOURCE_DIR}/regression/cases/meshes)

# Periodic conditions (step 8.3b): the table of the periodicity relations (both kinematics forms),
# the relations, T and T^T, the reduced operator and both diagonals against a serial assembly,
# the nodal weights, the zero-mean shift, Fbar and the volume-averaged stress (V8).
fepx_add_mpi_test(test_periodic SOURCES fe/test_periodic.f90 RANKS 1 2 4
                   ARGS ${CMAKE_CURRENT_SOURCE_DIR}/verification/cases/meshes)

# The isotropic viscoplastic start (D78): element matrices against v2's formulas, the exact
# affine solution of the three iterations on the periodic box.
fepx_add_mpi_test(test_isotropic_start SOURCES fe/test_isotropic_start.f90 RANKS 1 2 4
                   ARGS ${CMAKE_CURRENT_SOURCE_DIR}/verification/cases/meshes)
target_link_libraries(test_isotropic_start PRIVATE fepx_test_math_utils)

# One elastic increment (successive approximation), face forces and mid-node linearization.
fepx_add_mpi_test(test_increment SOURCES fe/test_increment.f90 RANKS 1 2 4
                   ARGS ${CMAKE_CURRENT_SOURCE_DIR}/regression/cases/meshes)

# Time stepping (step 5.2): rigid rotation of a pre-strained crystal (geometry update, lattice
# rotation by the continuum spin, elastic strain carried with the lattice).
fepx_add_mpi_test(test_rigid_rotation SOURCES fe/test_rigid_rotation.f90 RANKS 1 2 4
                   ARGS ${CMAKE_CURRENT_SOURCE_DIR}/regression/cases/meshes)

# Lamé verification (dev/theory/addendum_lame.md): quarter annulus, tet10 and hex20, two
# refinements, nu = 0.3, 0.49, 0.499 (tolerances: tests/verification/README.md).
fepx_add_mpi_test(test_lame SOURCES fe/test_lame.f90 RANKS 1 4
                   ARGS ${CMAKE_CURRENT_SOURCE_DIR}/verification/cases/meshes)

# .sim writer: every file of the step 4.3 subset read back and checked (identical on 1, 2 and 4
# ranks).
fepx_add_mpi_test(test_sim_writer SOURCES io/test_sim_writer.f90 RANKS 1 2 4
                   ARGS ${CMAKE_CURRENT_SOURCE_DIR}/regression/cases/meshes)

# Restart files (step 8.4): round trip of every record kind, the material state and a schedule
# step; error cases (header checks, records of another name, shape or type, a short file).
fepx_add_unit_test(test_restart_file SOURCES io/test_restart_file.f90)
# HDF5 files (step 9.1): only in a build with the HDF5 output.
if(FEPX_HAVE_HDF5)
  fepx_add_unit_test(test_hdf5_file SOURCES io/test_hdf5_file.f90)
endif()
add_executable(test_restart_file_errors io/test_restart_file_errors.f90)
target_link_libraries(test_restart_file_errors PRIVATE fepx_core)
set(restart_file_error_cases
    "not_restart|restart_open: errors.rst is not a FEPX restart file"
    "version|restart_open: errors.rst: restart format version 99, this FEPX reads version 2"
    "ranks|restart_open: errors.rst: written on 2 rank[(]s[)]. a restart needs the same number of ranks [(]this run: 3[)]"
    "rank|restart_open: errors.rst: the file of rank 1, read by rank 0"
    "missing|restart_open: cannot open restart file no_such_file.rst"
    "wrong_name|restart_get: errors.rst: record b expected, a found"
    "wrong_shape|restart_get: errors.rst: record a has shape [(]3[)], expected [(]4[)]"
    "wrong_type|restart_get: errors.rst: record a has another type"
    "end_of_file|restart_get: errors.rst: record a expected, end of file found"
    "state_size|restart_get_state: errors.rst: material state num_slip 2 in the restart file, 3 in this run")
foreach(error_case IN LISTS restart_file_error_cases)
  string(REPLACE "|" ";" error_case "${error_case}")
  list(GET error_case 0 case_name)
  list(GET error_case 1 case_message)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/restart_file_errors/${case_name})
  file(MAKE_DIRECTORY ${case_dir})
  add_test(NAME test_restart_file_error_${case_name}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=$<TARGET_FILE:test_restart_file_errors>;${case_name}"
                   -DEXPECTED_STATUS=1
                   "-DEXPECTED_STDOUT=Error  : ${case_message}"
                   "-DEXPECTED_STDERR=Error  : ${case_message}"
                   -DWORKING_DIR=${case_dir}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_restart_file_error_${case_name} PROPERTIES LABELS unit)
endforeach()

# .sim writer errors: every rank must end with the fatal error (MPI_Abort handler).
add_executable(test_sim_writer_errors io/test_sim_writer_errors.f90)
target_link_libraries(test_sim_writer_errors PRIVATE fepx_parallel)
set(sim_error_cases
    "dir_exists|sim_open: output directory err_dir_exists_np[0-9]+.sim already exists"
    "unknown_var|sim_open: internal error [(]please report it[)]: output variable strain_th is not written yet"
    "missing_input|sim_open: cannot read input file no_such_file.cfg"
    "node_shape|sim_write_nodes: internal error [(]please report it[)]: coo: expected 3 values per owned node"
    "elt_shape|sim_write_elts: internal error [(]please report it[)]: stress: expected 6 values per local element"
    "forces_shape|sim_write_forces: internal error [(]please report it[)]: expected a force and an area per face set"
    "slip_shape|sim_write_elts_slip: internal error [(]please report it[)]: crss: expected the values of every local element"
    "conv_shape|sim_write_convergence: internal error [(]please report it[)]: expected six norms and a CG count per iteration")
foreach(error_case IN LISTS sim_error_cases)
  string(REPLACE "|" ";" error_case "${error_case}")
  list(GET error_case 0 case_name)
  list(GET error_case 1 case_message)
  foreach(num_ranks 1 2)
    fepx_mpiexec_command(command ${num_ranks} $<TARGET_FILE:test_sim_writer_errors>)
    list(APPEND command ${case_name} ${CMAKE_CURRENT_SOURCE_DIR}/io/fixtures)
    add_test(NAME test_sim_writer_error_${case_name}_np${num_ranks}
             COMMAND ${CMAKE_COMMAND}
                     "-DCOMMAND=${command}"
                     -DEXPECTED_STATUS=nonzero
                     "-DEXPECTED_STDOUT=Error  : ${case_message}"
                     "-DEXPECTED_STDERR=Error  : ${case_message}"
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
    set_tests_properties(test_sim_writer_error_${case_name}_np${num_ranks} PROPERTIES
                         TIMEOUT ${fepx_timeout_60}
                         PROCESSORS ${num_ranks}
                         LABELS mpi)
  endforeach()
endforeach()

# Step schedule (step 5.1): times, strains, increments; and its errors, each a configuration
# written by the program to case.cfg in its own working directory (entries: name|message regex).
fepx_add_unit_test(test_step_schedule SOURCES fe/test_step_schedule.f90)
# Boundary-velocity iteration of triaxial load control, its pure parts (step 8.2b).
fepx_add_unit_test(test_bc_iteration SOURCES fe/test_bc_iteration.f90)
add_executable(test_step_schedule_errors fe/test_step_schedule_errors.f90)
target_link_libraries(test_step_schedule_errors PRIVATE fepx_core)
set(step_schedule_error_cases
    "rate_away|case.cfg: step 1, component 33: strain_rate 33 -1.0000E-03 moves away from target_strain 33 1.0000E-03 [(]strain at the start of the step: 0.0000E[+]00[)]"
    "reversal_rate_away|case.cfg: step 2, component 33: strain_rate 33 1.0000E-03 moves away from target_strain 33 1.0000E-03 [(]strain at the start of the step: 2.0000E-03[)]"
    "zero_change|case.cfg: step 2, component 33: target_strain 33 2.0000E-03 equals the strain at the start of the step"
    "time_not_after|case.cfg: step 2, component 33: target_time 1.5000E[+]00 s is not after the start of the step [(]2.0000E[+]00 s[)]"
    "dstrain_no_axis|case.cfg: step 1: dstrain needs a loading component"
    "dstrain_zero_rate|case.cfg: step 2, component 33: dstrain with a zero strain rate"
    "general_zero_vel|case.cfg: boundary_conditions general with a target on component 33: every set_bc velocity is zero"
    "dtime_min_dstrain|case.cfg: step 1, component 33: dtime_min 5.0000E-01 s is larger than the largest increment 1.0000E-01 s"
    "load_equal|case.cfg: steps 1 and 2, component 33: target_load 33 1.0000E[+]02 and 1.0000E[+]02 N are not more than twice load_tol [(]2 x 1.0000E-01 N[)] apart [(]a step ends within load_tol of its target[)]"
    "load_first_close|case.cfg: step 1, component 33: target_load 33 5.0000E-02 N is not more than load_tol 1.0000E-01 N from the initial value 0"
    "load_targets_close|case.cfg: steps 1 and 2, component 33: target_load 33 1.0000E[+]02 and 1.0015E[+]02 N are not more than twice load_tol [(]2 x 1.0000E-01 N[)] apart"
    "stress_targets_close|case.cfg: steps 1 and 2, component 33: target_stress 33 5.0000E[+]01 and 4.5000E[+]01 MPa are not more than twice stress_tol [(]2 x 3.0000E[+]00 MPa[)] apart"
    "load_away_previous|case.cfg: step 2, component 33: strain_rate 33 1.0000E-03 moves away from target_load 33 5.0000E[+]01 N [(]value at the start of the step: 1.0000E[+]02 N. a positive strain rate increases the load[)]"
    "rate_load_away|case.cfg: step 1, component 33: load_rate 33 -2.0000E[+]00 moves away from target_load 33 1.0000E[+]02 N [(]value at the start of the step: 0.0000E[+]00 N[)]"
    "rate_away_previous|case.cfg: step 2, component 33: load_rate 33 2.0000E[+]00 moves away from target_load 33 5.0000E[+]01 N [(]value at the start of the step: 1.0000E[+]02 N[)]"
    "rate_targets_close|case.cfg: steps 1 and 2, component 33: target_load 33 1.0000E[+]02 and 1.0015E[+]02 N are not more than twice load_tol"
    "rate_open_away|case.cfg: step 1, component 33: load_rate 33 -1.0000E[+]00 moves away from target_strain 33 1.0000E-03 [(]value at the start of the step: 0.0000E[+]00[)]"
    "max_incr_step|case.cfg: step 1, component 33: 2.0000E[+]01 increments of 1.0000E-01 s exceed max_incr 10"
    "max_incr_total|case.cfg: the steps need 1.1000E[+]01 increments, more than max_incr 10"
    "max_total_time|case.cfg: the steps end at 3.0000E[+]00 s, after max_total_time 2.5000E[+]00 s"
    "max_strain|case.cfg: step 2, component 33: the strain at the end of the step, -3.0000E-01, exceeds max_strain 2.0000E-01"
    "max_strain_time|case.cfg: step 1, component 33: the strain at the end of the step, 3.0000E-01, exceeds max_strain 2.0000E-01")
foreach(error_case IN LISTS step_schedule_error_cases)
  string(REPLACE "|" ";" error_case "${error_case}")
  list(GET error_case 0 case_name)
  list(GET error_case 1 case_message)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/step_schedule_errors/${case_name})
  file(MAKE_DIRECTORY ${case_dir})
  add_test(NAME test_step_schedule_error_${case_name}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=$<TARGET_FILE:test_step_schedule_errors>;${case_name}"
                   -DEXPECTED_STATUS=1
                   "-DEXPECTED_STDOUT=Error  : .*${case_message}"
                   "-DEXPECTED_STDERR=Error  : .*${case_message}"
                   -DWORKING_DIR=${case_dir}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_step_schedule_error_${case_name} PROPERTIES LABELS unit)
endforeach()

# Boundary-condition errors: every rank must end with the fatal error (MPI_Abort handler).
add_executable(test_fe_errors fe/test_fe_errors.f90)
target_link_libraries(test_fe_errors PRIVATE fepx_parallel)
set(fe_error_cases
    "unknown_nset|bc_resolve: node set x9 is not in the mesh"
    "conflict|bc_resolve: conflicting velocities prescribed on a node of node set face \\(component x\\)"
    "bad_axis|bc_preset_specs: internal error [(]please report it[)]: the loading axis must be 1, 2 or 3 \\(x, y or z\\)"
    "bad_component|bc_resolve: internal error [(]please report it[)]: velocity component of node set apex must be 1, 2 or 3 \\(x, y or z\\)"
    "mpc_unknown_nset|constraint_table_mpc: line 1: node set x9 of a multi-point constraint is not in the mesh"
    "mpc_two_lines|constraint_table_mpc: node 3 component x is a secondary of two different multi-point constraints \\(lines 1 and 2\\)"
    "mpc_cycle|constraint_table_mpc: the multi-point constraints of line [12] form a cycle"
    "mpc_conflict|constraint_propagate_bcs: conflicting velocity conditions in a multi-point constraint group")
foreach(error_case IN LISTS fe_error_cases)
  string(REPLACE "|" ";" error_case "${error_case}")
  list(GET error_case 0 case_name)
  list(GET error_case 1 case_message)
  foreach(num_ranks 1 2)
    fepx_mpiexec_command(command ${num_ranks} $<TARGET_FILE:test_fe_errors>)
    list(APPEND command ${case_name} ${CMAKE_CURRENT_SOURCE_DIR}/io/fixtures)
    add_test(NAME test_fe_error_${case_name}_np${num_ranks}
             COMMAND ${CMAKE_COMMAND}
                     "-DCOMMAND=${command}"
                     -DEXPECTED_STATUS=nonzero
                     "-DEXPECTED_STDOUT=Error  : ${case_message}"
                     "-DEXPECTED_STDERR=Error  : ${case_message}"
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
    set_tests_properties(test_fe_error_${case_name}_np${num_ranks} PROPERTIES
                         TIMEOUT ${fepx_timeout_60}
                         PROCESSORS ${num_ranks}
                         LABELS mpi)
  endforeach()
endforeach()

# fepx restart errors (step 8.4): part 1 of verification/cases/restart_monotonic on 2 ranks,
# then part 2 on another rank count, with a changed material, without the
# output directory, without restart on, with other output formats (each an error with its
# message), with output_format ascii given (the default: no change), and with a changed
# solver option (a warning; tests/app/restart_errors.cmake).
add_test(NAME test_fepx_restart_errors
         COMMAND ${CMAKE_COMMAND}
                 -DEXECUTABLE=$<TARGET_FILE:fepx>
                 -DCASE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/verification/cases/restart_monotonic
                 -DMESH=${CMAKE_CURRENT_SOURCE_DIR}/regression/cases/meshes/sx2_part2.msh
                 -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/fepx_restart_errors
                 -DMPIEXEC=${MPIEXEC_EXECUTABLE}
                 -DMPIEXEC_NUMPROC=${MPIEXEC_NUMPROC_FLAG}
                 "-DMPIEXEC_EXTRA=${FEPX_MPIEXEC_EXTRA_FLAGS_TEXT}"
                 -DHDF5=${FEPX_HAVE_HDF5}
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/app/restart_errors.cmake)
set_tests_properties(test_fepx_restart_errors PROPERTIES
                     PROCESSORS 2
                     TIMEOUT ${fepx_timeout_600}
                     LABELS "mpi;fepx")

# fepx executable errors: each case is a configuration (app/cases/<case>.cfg) run on the sx2 mesh
# in its own working directory, on 1 and 2 ranks; every rank must end with the fatal error.
set(app_error_cases
    "phase_count|fepx: the mesh has 1 phase[(]s[)], the configuration 2 [(]number_of_phases[)]"
    "rate_away|cfg: simulation.cfg: step 1, component 33: strain_rate 33 -1.0000E-03 moves away from target_strain 33 1.0000E-04"
    "no_target|cfg: simulation.cfg: step 1, component 33: no target [(]target_strain 33, target_load 33, target_stress 33, target_time or target_duration[)]"
    "load_rate_away|cfg: simulation.cfg: step 1, component 33: strain_rate 33 -1.0000E-03 moves away from target_load 33 1.0000E[+]02 N [(]value at the start of the step: 0.0000E[+]00 N. a positive strain rate increases the load[)]"
    "stress_rate_away_run|cfg: simulation.cfg: step 2, component 33: strain_rate 33 1.0000E-03 moves away from target_stress 33 0.0000E[+]00 MPa [(]value at the start of the step: 1.[0-9]+E[+]01 MPa"
    "stress_reached_start|cfg: simulation.cfg: step 2, component 33: target_stress 33 1.5000E[+]01 MPa is already reached at the start of the step [(]1.[0-9]+E[+]01 MPa[)]"
    "scys_vertex_rank|fepx: initialization isotropic_viscoplastic: the active slip systems of phase 1 span fewer than 5 deviatoric dimensions [(]at least five independent systems are needed[)]"
    "opt_both_forms|fepx: simulation.opt: give [$]ElsetCrss or [$]ElementCrss, not both"
    "opt_count|fepx: simulation.opt: Crss of elset 1 [(]phase 1[)]: 2 value[(]s[)]. give one value, one per active slip family or one per slip system"
    "opt_tau_s_low|fepx: simulation.opt: element [0-9]+: the saturation strength must exceed the initial strength"
    "km_opt_both_forms|fepx: simulation.opt: give [$]ElsetRho0 or [$]ElementRho0, not both"
    "km_opt_count|fepx: simulation.opt: Rho0 of elset 1 [(]phase 1[)]: 2 value[(]s[)]. give one value"
    "km_opt_zero|fepx: simulation.opt: elset 1: the initial dislocation density Rho0 must be positive"
    "irr_opt_both_forms|fepx: simulation.opt: give [$]ElsetIrrLoopDensity or [$]ElementIrrLoopDensity, not both"
    "irr_opt_negative|fepx: simulation.opt: elset 1: the loop density IrrLoopDensity must not be negative"
    "irr_opt_diameter_zero|fepx: simulation.opt: elset 1: the loop diameter IrrLoopDiameter must be positive"
    "subincr_depth|fepx: increment 1: 2175 point solve[(]s[)] did not converge at the converged velocities after 2 subdivisions [(]sub-increment dt/2.2 =  1.22681E-01 s[)]. refine the time stepping [(]number_of_incrs, dtime or dstrain[)] or raise subincrement_max_depth"
    "subincr_successive|fepx: 2 successive increments needed sub-incrementation [(]increment 3. the global iteration did not converge in global_max_iters iterations[)]. refine the time stepping [(]number_of_incrs, dtime or dstrain[)] or raise subincrement_max_successive")
foreach(error_case IN LISTS app_error_cases)
  string(REPLACE "|" ";" error_case "${error_case}")
  list(GET error_case 0 case_name)
  list(GET error_case 1 case_message)
  foreach(num_ranks 1 2)
    set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/app_errors/${case_name}_np${num_ranks})
    file(MAKE_DIRECTORY ${case_dir})
    configure_file(app/cases/${case_name}.cfg ${case_dir}/simulation.cfg COPYONLY)
    fepx_mesh_for_ranks(case_mesh regression/cases/meshes/sx2.msh ${num_ranks})
    configure_file(${case_mesh} ${case_dir}/simulation.msh COPYONLY)
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/app/cases/${case_name}.opt)
      configure_file(app/cases/${case_name}.opt ${case_dir}/simulation.opt COPYONLY)
    endif()
    fepx_mpiexec_command(command ${num_ranks} $<TARGET_FILE:fepx>)
    add_test(NAME test_fepx_error_${case_name}_np${num_ranks}
             COMMAND ${CMAKE_COMMAND}
                     "-DCOMMAND=${command}"
                     -DEXPECTED_STATUS=nonzero
                     "-DEXPECTED_STDOUT=Error  : ${case_message}"
                     "-DEXPECTED_STDERR=Error  : ${case_message}"
                     -DWORKING_DIR=${case_dir}
                     -DREMOVE_BEFORE=${case_dir}/simulation.sim
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
    set_tests_properties(test_fepx_error_${case_name}_np${num_ranks} PROPERTIES
                         TIMEOUT ${fepx_timeout_60}
                         PROCESSORS ${num_ranks}
                         LABELS mpi)
  endforeach()
endforeach()

# fepx executable errors with periodic conditions (step 8.3b, addendum_periodic_mpc.md §5 V7):
# each case on its own mesh (case|mesh|message), on 1 and 2 ranks. periodic_free_rotation also
# checks the three free-rotation warnings (Q-P2) before the divergence (v2's behaviour, P2).
set(periodic_error_cases
    "periodic_no_relations|regression/cases/meshes/sx2.msh|periodic: boundary_conditions periodic needs a mesh with periodicity relations"
    "periodic_axis_x|verification/cases/meshes/periodic_cube_x.msh|periodic: the mesh is not periodic along y"
    "periodic_bc_secondary|verification/cases/meshes/periodic_box_2x1x1.msh|constraint_check_bcs: a velocity condition is on a periodic secondary dof"
    "periodic_free_rotation|verification/cases/meshes/periodic_cube.msh|fepx: the global iteration did not converge in global_max_iters iterations")
foreach(error_case IN LISTS periodic_error_cases)
  string(REPLACE "|" ";" error_case "${error_case}")
  list(GET error_case 0 case_name)
  list(GET error_case 1 case_mesh)
  list(GET error_case 2 case_message)
  set(prefix "")
  if(case_name STREQUAL "periodic_free_rotation")
    set(prefix "Periodic conditions: Fbar-dot 12 and 21 are both free.*Fbar-dot 13 and 31.*Fbar-dot 23 and 32.*")
  endif()
  foreach(num_ranks 1 2)
    set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/app_errors/${case_name}_np${num_ranks})
    file(MAKE_DIRECTORY ${case_dir})
    configure_file(app/cases/${case_name}.cfg ${case_dir}/simulation.cfg COPYONLY)
    fepx_mesh_for_ranks(case_mesh_ranks ${case_mesh} ${num_ranks})
    configure_file(${case_mesh_ranks} ${case_dir}/simulation.msh COPYONLY)
    fepx_mpiexec_command(command ${num_ranks} $<TARGET_FILE:fepx>)
    add_test(NAME test_fepx_error_${case_name}_np${num_ranks}
             COMMAND ${CMAKE_COMMAND}
                     "-DCOMMAND=${command}"
                     -DEXPECTED_STATUS=nonzero
                     "-DEXPECTED_STDOUT=${prefix}Error  : ${case_message}"
                     "-DEXPECTED_STDERR=Error  : ${case_message}"
                     -DWORKING_DIR=${case_dir}
                     -DREMOVE_BEFORE=${case_dir}/simulation.sim
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
    set_tests_properties(test_fepx_error_${case_name}_np${num_ranks} PROPERTIES
                         TIMEOUT ${fepx_timeout_60}
                         PROCESSORS ${num_ranks}
                         LABELS mpi)
  endforeach()
endforeach()

# fepx executable run of a crystal-plasticity phase (step 6.2, D49): the initialization elastic
# step (0.01 * 3 tau_0 / max row sum of [C] / rate = 0.03 * 100 / 323400 / 1e-3 s) and a first
# increment shortened by it; the run completes.
foreach(num_ranks 1 2)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/app_runs/crystal_plasticity_np${num_ranks})
  file(MAKE_DIRECTORY ${case_dir})
  configure_file(app/cases/crystal_plasticity_run.cfg ${case_dir}/simulation.cfg COPYONLY)
  fepx_mesh_for_ranks(case_mesh regression/cases/meshes/sx2.msh ${num_ranks})
  configure_file(${case_mesh} ${case_dir}/simulation.msh COPYONLY)
  fepx_mpiexec_command(command ${num_ranks} $<TARGET_FILE:fepx>)
  add_test(NAME test_fepx_run_crystal_plasticity_np${num_ranks}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=${command}"
                   -DEXPECTED_STATUS=0
                   "-DEXPECTED_STDOUT=Elastic step of 9.2764E-03 s [(]step 1[)].*Running step 1[.][.][.].*Increment 1: t = [0-9.]+E[-+][0-9]+ s, dt = 9.9072E-01 s.*Simulation completed successfully"
                   -DWORKING_DIR=${case_dir}
                   -DREMOVE_BEFORE=${case_dir}/simulation.sim
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_fepx_run_crystal_plasticity_np${num_ranks} PROPERTIES
                       TIMEOUT ${fepx_timeout_120}
                       PROCESSORS ${num_ranks}
                       LABELS mpi)
endforeach()

# fepx executable run with simulation.opt strengths on a two-phase mesh (step 7.1, D59): phase 1
# (crystal plasticity, precipitate_strengthening on) and phase 2 (elastic, via simulation.phase);
# both warnings are logged (precipitates ignored for strengths from the file; values given for an
# elastic phase ignored) and the run completes.
foreach(num_ranks 1 2)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/app_runs/opt_warnings_np${num_ranks})
  file(MAKE_DIRECTORY ${case_dir})
  configure_file(app/cases/opt_warnings_run.cfg ${case_dir}/simulation.cfg COPYONLY)
  configure_file(app/cases/opt_warnings_run.opt ${case_dir}/simulation.opt COPYONLY)
  configure_file(app/cases/opt_warnings_run.phase ${case_dir}/simulation.phase COPYONLY)
  fepx_mesh_for_ranks(case_mesh regression/cases/meshes/sx2.msh ${num_ranks})
  configure_file(${case_mesh} ${case_dir}/simulation.msh COPYONLY)
  fepx_mpiexec_command(command ${num_ranks} $<TARGET_FILE:fepx>)
  add_test(NAME test_fepx_run_opt_warnings_np${num_ranks}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=${command}"
                   -DEXPECTED_STATUS=0
                   "-DEXPECTED_STDOUT=simulation.opt gives plasticity parameters for phase 2, which is elastic: ignored.*Phase 1: .*Initial strengths from simulation.opt are used as given. precipitate_strengthening is ignored for them.*Simulation completed successfully"
                   -DWORKING_DIR=${case_dir}
                   -DREMOVE_BEFORE=${case_dir}/simulation.sim
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_fepx_run_opt_warnings_np${num_ranks} PROPERTIES
                       TIMEOUT ${fepx_timeout_120}
                       PROCESSORS ${num_ranks}
                       LABELS mpi)
endforeach()

# fepx executable run with simulation.opt values on a two-phase mesh with Kocks-Mecking hardening
# (step 7.2, D60): phase 1 Kocks-Mecking (only $ElsetRho0 is read), phase 2 saturation_isotropic
# (Rho0 not used); the file gives $ElsetCrss and $ElsetRho0 for both; both warnings are logged,
# with Hill's mu_0 of phase 1, and the run completes.
foreach(num_ranks 1 2)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/app_runs/km_opt_warnings_np${num_ranks})
  file(MAKE_DIRECTORY ${case_dir})
  configure_file(app/cases/km_opt_warnings_run.cfg ${case_dir}/simulation.cfg COPYONLY)
  configure_file(app/cases/km_opt_warnings_run.opt ${case_dir}/simulation.opt COPYONLY)
  configure_file(app/cases/km_opt_warnings_run.phase ${case_dir}/simulation.phase COPYONLY)
  fepx_mesh_for_ranks(case_mesh regression/cases/meshes/sx2.msh ${num_ranks})
  configure_file(${case_mesh} ${case_dir}/simulation.msh COPYONLY)
  fepx_mpiexec_command(command ${num_ranks} $<TARGET_FILE:fepx>)
  add_test(NAME test_fepx_run_km_opt_warnings_np${num_ranks}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=${command}"
                   -DEXPECTED_STATUS=0
                   "-DEXPECTED_STDOUT=simulation.opt gives Crss, CrssSat or HardCoeff values for phase 1, whose hardening is kocks_mecking [(]only Rho0 and the irradiation loop sections are read[)]: ignored.*simulation.opt gives Rho0 values for phase 2, whose hardening is not kocks_mecking: ignored.*Phase 1: .*Kocks-Mecking shear modulus mu_0 [(]Hill[)]: 5.479[0-9]E[+]04 MPa.*Simulation completed successfully"
                   -DWORKING_DIR=${case_dir}
                   -DREMOVE_BEFORE=${case_dir}/simulation.sim
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_fepx_run_km_opt_warnings_np${num_ranks} PROPERTIES
                       TIMEOUT ${fepx_timeout_120}
                       PROCESSORS ${num_ranks}
                       LABELS mpi)
endforeach()

# fepx executable run with sub-incrementation (Q35, D63): point_max_iters_state 3 makes the
# point solves of increments 2-4 fail at full length; each is solved as halves, subdivided where
# rejected (increment 2: 5 sub-increments in 9 attempts down to dt/16, 3: 9 in 17, 4: 4 in 7, the
# same on 1 and 2 ranks); the convergence file has the attempt columns and the run completes.
foreach(num_ranks 1 2)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/app_runs/subincr_np${num_ranks})
  file(MAKE_DIRECTORY ${case_dir})
  configure_file(app/cases/subincr_run.cfg ${case_dir}/simulation.cfg COPYONLY)
  fepx_mesh_for_ranks(case_mesh regression/cases/meshes/sx2.msh ${num_ranks})
  configure_file(${case_mesh} ${case_dir}/simulation.msh COPYONLY)
  fepx_mpiexec_command(command ${num_ranks} $<TARGET_FILE:fepx>)
  add_test(NAME test_fepx_run_subincr_np${num_ranks}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=${command}"
                   -DEXPECTED_STATUS=0
                   "-DEXPECTED_STDOUT=Increment 2: .*Attempt over 1.0000E[+]00 s rejected [(]2175 point solve[(]s[)] did not converge at the converged velocities[)]. subdividing [(]depth 1[)].*5 sub-increments [(]9 attempts[)], smallest 6.2500E-02 s.*Increment 3: .*9 sub-increments [(]17 attempts[)], smallest 6.2500E-02 s.*Increment 4: .*4 sub-increments [(]7 attempts[)], smallest 2.5000E-01 s.*Simulation completed successfully"
                   -DWORKING_DIR=${case_dir}
                   -DREMOVE_BEFORE=${case_dir}/simulation.sim
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_fepx_run_subincr_np${num_ranks} PROPERTIES
                       TIMEOUT ${fepx_timeout_120}
                       PROCESSORS ${num_ranks}
                       LABELS mpi)
endforeach()

# fepx executable run of two steps (step 5.2): 4 increments, then 4 with a shortened last one
# (dtime 0.03 s in a 0.1 s step); the last increment ends at t = 0.3 s and the run completes.
foreach(num_ranks 1 2)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/app_runs/multi_step_np${num_ranks})
  file(MAKE_DIRECTORY ${case_dir})
  configure_file(app/cases/multi_step_run.cfg ${case_dir}/simulation.cfg COPYONLY)
  fepx_mesh_for_ranks(case_mesh regression/cases/meshes/sx2.msh ${num_ranks})
  configure_file(${case_mesh} ${case_dir}/simulation.msh COPYONLY)
  fepx_mpiexec_command(command ${num_ranks} $<TARGET_FILE:fepx>)
  add_test(NAME test_fepx_run_multi_step_np${num_ranks}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=${command}"
                   -DEXPECTED_STATUS=0
                   "-DEXPECTED_STDOUT=Step 2: t = 2.0000E-01 -> 3.0000E-01 s, 4 increments.*Running step 2[.][.][.].*Increment 8: t = [0-9.]+E[-+][0-9]+ s, dt = 1.0000E-02 s.*Simulation completed successfully"
                   -DWORKING_DIR=${case_dir}
                   -DREMOVE_BEFORE=${case_dir}/simulation.sim
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_fepx_run_multi_step_np${num_ranks} PROPERTIES
                       TIMEOUT ${fepx_timeout_120}
                       PROCESSORS ${num_ranks}
                       LABELS mpi)
endforeach()

# Command line of the fepx program (step 9.7, Q64): --version and --help are handled before MPI
# starts (run without the launcher) and end with status 0; any other argument is an error
# (status 1, the message on both units).
add_test(NAME test_fepx_version
         COMMAND ${CMAKE_COMMAND}
                 "-DCOMMAND=$<TARGET_FILE:fepx>;--version"
                 -DEXPECTED_STATUS=0
                 "-DEXPECTED_STDOUT=^[0-9]+[.][0-9]+[.][0-9]+(-dev|[+][0-9]+)?( [(]git [0-9a-f]+(-dirty)?, [0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9][)])?"
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
add_test(NAME test_fepx_help
         COMMAND ${CMAKE_COMMAND}
                 "-DCOMMAND=$<TARGET_FILE:fepx>;--help"
                 -DEXPECTED_STATUS=0
                 "-DEXPECTED_STDOUT=^=+    F   E   P   X    =+.*Info   : Version .*Usage: mpirun -np <ranks> fepx .*fepx --version"
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
add_test(NAME test_fepx_bad_option
         COMMAND ${CMAKE_COMMAND}
                 "-DCOMMAND=$<TARGET_FILE:fepx>;--bogus"
                 -DEXPECTED_STATUS=1
                 "-DEXPECTED_STDOUT=^Error  : unknown option \"--bogus\" [(]see fepx --help[)]"
                 "-DEXPECTED_STDERR=^Error  : unknown option \"--bogus\" [(]see fepx --help[)]"
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
# The manual page (step 9.8, D101) describes every option of the help text and formats cleanly.
find_program(FEPX_GROFF_PROGRAM groff)
if(NOT FEPX_GROFF_PROGRAM)
  set(FEPX_GROFF_PROGRAM "")
endif()
add_test(NAME test_fepx_man_page
         COMMAND ${CMAKE_COMMAND}
                 -DEXECUTABLE=$<TARGET_FILE:fepx>
                 -DMAN_PAGE=${FEPX_ROOT_DIR}/doc/man/fepx.1
                 -DGROFF=${FEPX_GROFF_PROGRAM}
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/app/man_page.cmake)
set_tests_properties(test_fepx_version test_fepx_help test_fepx_bad_option test_fepx_man_page
                     PROPERTIES TIMEOUT ${fepx_timeout_60})

# The number of ranks and the mesh's partitions (D108): on more than one rank they must agree (a
# mesh of another number of partitions and an unpartitioned mesh are errors); 1 rank takes any
# mesh, with a warning when it is partitioned.
foreach(partition_case
        "mismatch|regression/cases/meshes/poly_part4.msh|2|nonzero|Error  : mesh_distribute: The mesh has 4 partition[(]s[)], but 2 ranks were started: partition it for 2 ranks"
        "unpartitioned|regression/cases/meshes/sx2.msh|2|nonzero|Error  : mesh_distribute: The mesh is not partitioned, but 2 ranks were started"
        "one_rank|regression/cases/meshes/poly_part4.msh|1|0|Warning:     > The mesh has 4 partitions, but 1 rank was started: the rank takes")
  string(REPLACE "|" ";" partition_case "${partition_case}")
  list(GET partition_case 0 case_name)
  list(GET partition_case 1 case_mesh)
  list(GET partition_case 2 case_ranks)
  list(GET partition_case 3 case_status)
  list(GET partition_case 4 case_message)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/app_runs/partition_${case_name})
  file(MAKE_DIRECTORY ${case_dir})
  configure_file(app/cases/crystal_plasticity_run.cfg ${case_dir}/simulation.cfg COPYONLY)
  configure_file(${case_mesh} ${case_dir}/simulation.msh COPYONLY)
  fepx_mpiexec_command(command ${case_ranks} $<TARGET_FILE:fepx>)
  add_test(NAME test_fepx_partition_${case_name}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=${command}"
                   -DEXPECTED_STATUS=${case_status}
                   "-DEXPECTED_STDOUT=${case_message}"
                   -DWORKING_DIR=${case_dir}
                   -DREMOVE_BEFORE=${case_dir}/simulation.sim
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_fepx_partition_${case_name} PROPERTIES
                       TIMEOUT ${fepx_timeout_60}
                       PROCESSORS ${case_ranks}
                       LABELS mpi)
endforeach()

# An input error that every rank detects (an unknown keyword), on 2 ranks (step 9.7, Q62): the
# message once on standard output (rank 0's; another rank's would be tagged "[Rank k]"), and on
# standard error.
set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/app_errors/error_once_np2)
file(MAKE_DIRECTORY ${case_dir})
file(READ app/cases/crystal_plasticity_run.cfg error_once_cfg)
file(WRITE ${case_dir}/simulation.cfg "${error_once_cfg}bogus_keyword 1\n")
fepx_mesh_for_ranks(case_mesh regression/cases/meshes/sx2.msh 2)
configure_file(${case_mesh} ${case_dir}/simulation.msh COPYONLY)
fepx_mpiexec_command(command 2 $<TARGET_FILE:fepx>)
add_test(NAME test_fepx_error_once_np2
         COMMAND ${CMAKE_COMMAND}
                 "-DCOMMAND=${command}"
                 -DEXPECTED_STATUS=nonzero
                 "-DEXPECTED_STDOUT_ONCE=Error  : cfg: simulation.cfg:[0-9]+: unknown keyword"
                 "-DEXPECTED_STDERR=Error  : cfg: simulation.cfg:[0-9]+: unknown keyword"
                 -DWORKING_DIR=${case_dir}
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
set_tests_properties(test_fepx_error_once_np2 PROPERTIES
                     TIMEOUT ${fepx_timeout_60}
                     PROCESSORS 2
                     LABELS mpi)

# fepx executable run with general conditions and a compressive strain target (step 5.1): the
# step lasts |de| L / max |v| = 1e-4 * 1 / 1e-3 = 0.1 s at strain rate -1e-3/s (set_bc velocities
# reversed); one increment, completed.
foreach(num_ranks 1 2)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/app_runs/general_strain_target_np${num_ranks})
  file(MAKE_DIRECTORY ${case_dir})
  configure_file(app/cases/general_strain_target_run.cfg ${case_dir}/simulation.cfg COPYONLY)
  fepx_mesh_for_ranks(case_mesh regression/cases/meshes/sx2.msh ${num_ranks})
  configure_file(${case_mesh} ${case_dir}/simulation.msh COPYONLY)
  fepx_mpiexec_command(command ${num_ranks} $<TARGET_FILE:fepx>)
  add_test(NAME test_fepx_run_general_strain_target_np${num_ranks}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=${command}"
                   -DEXPECTED_STATUS=0
                   "-DEXPECTED_STDOUT=Step 1: t = 0.0000E[+]00 -> 1.0000E-01 s, 1 increment.*Strain 33: 0.0000E[+]00 -> -1.0000E-04 at -1.0000E-03/s.*Simulation completed successfully"
                   -DWORKING_DIR=${case_dir}
                   -DREMOVE_BEFORE=${case_dir}/simulation.sim
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_fepx_run_general_strain_target_np${num_ranks} PROPERTIES
                       TIMEOUT ${fepx_timeout_60}
                       PROCESSORS ${num_ranks}
                       LABELS mpi)
endforeach()

# fepx executable runs with load and stress targets (step 8.1, addendum_load_targets.md §3):
#   load_unload_run   elastic: a stress target, unloading to zero load, then a strain target
#                     planned during the run from the strain reached;
#   necking           the v2-fixed C54 case started as v2 (isotropic_viscoplastic): the load peaks at
#                     increment 20 and falls at 21 (33 and 34 on the sx2 mesh of Neper 713); the run stops there with the
#                     error and writes step 1 from that increment (V6).
foreach(num_ranks 1 2)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/app_runs/load_unload_np${num_ranks})
  file(MAKE_DIRECTORY ${case_dir})
  configure_file(app/cases/load_unload_run.cfg ${case_dir}/simulation.cfg COPYONLY)
  fepx_mesh_for_ranks(case_mesh regression/cases/meshes/sx2.msh ${num_ranks})
  configure_file(${case_mesh} ${case_dir}/simulation.msh COPYONLY)
  fepx_mpiexec_command(command ${num_ranks} $<TARGET_FILE:fepx>)
  add_test(NAME test_fepx_run_load_unload_np${num_ranks}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=${command}"
                   -DEXPECTED_STATUS=0
                   "-DEXPECTED_STDOUT=Step 1: target_stress 33 4.0000E[+]01 MPa from t = 0.0000E[+]00 s.*Step 1: target 4.000000E[+]01, reached 4.0[0-9]+E[+]01.*Step 2: target_load 33 0.0000E[+]00 N from t = [0-9.]+E-01 s.*Step 2: target 0.000000E[+]00, reached -[0-9.]+E-0[0-9].*Step 3: t = [0-9.]+E-01 -> [0-9.]+E-01 s, 2 increments.*Strain 33: -?[0-9.]+E-0[6-9] -> 1.0000E-04 at 1.0000E-03/s.*Simulation completed successfully"
                   -DWORKING_DIR=${case_dir}
                   -DREMOVE_BEFORE=${case_dir}/simulation.sim
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_fepx_run_load_unload_np${num_ranks} PROPERTIES
                       TIMEOUT ${fepx_timeout_60}
                       PROCESSORS ${num_ranks}
                       LABELS mpi)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/app_runs/necking_np${num_ranks})
  file(MAKE_DIRECTORY ${case_dir})
  configure_file(app/cases/necking_run.cfg ${case_dir}/simulation.cfg COPYONLY)
  fepx_mesh_for_ranks(case_mesh regression/cases/meshes/sx2.msh ${num_ranks})
  configure_file(${case_mesh} ${case_dir}/simulation.msh COPYONLY)
  add_test(NAME test_fepx_run_necking_np${num_ranks}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=${command}"
                   -DEXPECTED_STATUS=nonzero
                   "-DEXPECTED_STDOUT=Error  : fepx: stopped in step 1 at increment 21 [(]t =  2.0997E[+]01 s[)]: necking [(]check_necking on[)]: the load on face z1 fell from  6.1291[0-9]+E[+]02 N to  6.1236[0-9]+E[+]02 N"
                   -DEXPECTED_FILE=simulation.sim/results/elts/stress/stress.step1
                   -DWORKING_DIR=${case_dir}
                   -DREMOVE_BEFORE=${case_dir}/simulation.sim
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_fepx_run_necking_np${num_ranks} PROPERTIES
                       TIMEOUT ${fepx_timeout_120}
                       PROCESSORS ${num_ranks}
                       LABELS mpi)
endforeach()

# fepx executable runs of triaxial load control (step 8.2b, addendum_triaxial_load_control.md §3.4,
# §5 V6):
#   subincr_triaxial_run  hcp crystal, control 33 to 0.4 % with stress 11 in time to -20 MPa,
#                         subincrementation on with point_max_iters_state 3 (as subincr_run):
#                         solves of the boundary-velocity iteration, perturbed ones included,
#                         reject attempts, which are subdivided; the run completes;
#   bc_iter_max           the isotropic elastic cube of verification_fepx_iso_triaxial_csr with
#                         max_bc_iter 1: the first increment needs a correction, the run stops.
foreach(num_ranks 1 2)
  fepx_mpiexec_command(command ${num_ranks} $<TARGET_FILE:fepx>)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/app_runs/subincr_triaxial_np${num_ranks})
  file(MAKE_DIRECTORY ${case_dir})
  configure_file(app/cases/subincr_triaxial_run.cfg ${case_dir}/simulation.cfg COPYONLY)
  fepx_mesh_for_ranks(case_mesh regression/cases/meshes/sx2.msh ${num_ranks})
  configure_file(${case_mesh} ${case_dir}/simulation.msh COPYONLY)
  add_test(NAME test_fepx_run_subincr_triaxial_np${num_ranks}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=${command}"
                   -DEXPECTED_STATUS=0
                   "-DEXPECTED_STDOUT=Load iteration 1: .*Increment 2: .*Attempt over 1.0000E[+]00 s rejected .*subdividing [(]depth 1[)].*15 sub-increments [(]29 attempts[)].*Increment 3: .*38 sub-increments [(]75 attempts[)].*Increment 4: .*16 sub-increments [(]31 attempts[)].*Simulation completed successfully"
                   -DWORKING_DIR=${case_dir}
                   -DREMOVE_BEFORE=${case_dir}/simulation.sim
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_fepx_run_subincr_triaxial_np${num_ranks} PROPERTIES
                       TIMEOUT ${fepx_timeout_300}
                       PROCESSORS ${num_ranks}
                       LABELS mpi)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/app_runs/bc_iter_max_np${num_ranks})
  file(MAKE_DIRECTORY ${case_dir})
  configure_file(app/cases/bc_iter_max.cfg ${case_dir}/simulation.cfg COPYONLY)
  fepx_mesh_for_ranks(case_mesh verification/cases/meshes/sx2_id.msh ${num_ranks})
  configure_file(${case_mesh} ${case_dir}/simulation.msh COPYONLY)
  add_test(NAME test_fepx_run_bc_iter_max_np${num_ranks}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=${command}"
                   -DEXPECTED_STATUS=nonzero
                   "-DEXPECTED_STDOUT=Error  : fepx: step 1, increment 1: the boundary-velocity iteration did not converge in max_bc_iter 1 iterations [(]component 22: [ -][0-9.]+E[-+][0-9]+, to reach  0.000000E[+]00 N. component 33:  4.99[0-9]+E[+]01, to reach  5.01[0-9]+E[+]01 N[)]. raise max_bc_iter or the tolerances [(]load_tol, stress_tol[)], or shorten the increments"
                   -DWORKING_DIR=${case_dir}
                   -DREMOVE_BEFORE=${case_dir}/simulation.sim
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_fepx_run_bc_iter_max_np${num_ranks} PROPERTIES
                       TIMEOUT ${fepx_timeout_120}
                       PROCESSORS ${num_ranks}
                       LABELS mpi)
endforeach()

# fepx executable run of iteration (B) beyond max_bc_iter (step 8.2c, addendum §3.4, §5 V6): the
# creep case of verification_fepx_fcc_sx_creep_dwell with max_bc_iter 1; the first increment's
# scaling loop needs more iterations, the run stops with the error.
foreach(num_ranks 1 2)
  fepx_mpiexec_command(command ${num_ranks} $<TARGET_FILE:fepx>)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/app_runs/bc_rate_iter_max_np${num_ranks})
  file(MAKE_DIRECTORY ${case_dir})
  configure_file(app/cases/bc_rate_iter_max.cfg ${case_dir}/simulation.cfg COPYONLY)
  fepx_mesh_for_ranks(case_mesh regression/cases/meshes/sx2.msh ${num_ranks})
  configure_file(${case_mesh} ${case_dir}/simulation.msh COPYONLY)
  add_test(NAME test_fepx_run_bc_rate_iter_max_np${num_ranks}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=${command}"
                   -DEXPECTED_STATUS=nonzero
                   "-DEXPECTED_STDOUT=Error  : fepx: step 1, increment 1: the boundary-velocity iteration [(]scaling loop[)] did not converge in max_bc_iter 1 iterations [(]component 33:  6.[0-9]+E[+]01, to reach  5.500000E[+]01 N[)]. raise max_bc_iter"
                   -DWORKING_DIR=${case_dir}
                   -DREMOVE_BEFORE=${case_dir}/simulation.sim
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_fepx_run_bc_rate_iter_max_np${num_ranks} PROPERTIES
                       TIMEOUT ${fepx_timeout_120}
                       PROCESSORS ${num_ranks}
                       LABELS mpi)
endforeach()

# fepx executable runs landing on load targets (step 8.2a, addendum_triaxial_load_control.md §3.6,
# §5 V6-V7):
#   landing_run          elastic-plastic fcc crystal to 400 N near yield and back to 0 N with
#                        dtime_factor 2: increment 3 lands in 5 attempts (false position on the
#                        bending response), the unloading in 1;
#   landing_max_bc_iter  the same with max_bc_iter 3: the run stops with the error after writing
#                        the last accepted increment as step 1.
foreach(num_ranks 1 2)
  fepx_mpiexec_command(command ${num_ranks} $<TARGET_FILE:fepx>)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/app_runs/landing_np${num_ranks})
  file(MAKE_DIRECTORY ${case_dir})
  configure_file(app/cases/landing_run.cfg ${case_dir}/simulation.cfg COPYONLY)
  fepx_mesh_for_ranks(case_mesh regression/cases/meshes/sx2.msh ${num_ranks})
  configure_file(${case_mesh} ${case_dir}/simulation.msh COPYONLY)
  add_test(NAME test_fepx_run_landing_np${num_ranks}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=${command}"
                   -DEXPECTED_STATUS=0
                   "-DEXPECTED_STDOUT=Increment 3: .*Landed at 4.000[0-9]+E[+]02 N in 5 attempts.*Step 1: target 4.000000E[+]02, reached 4.000[0-9]+E[+]02.*Increment 6: .*Landed at -?[0-9.]+E-0[0-9] N in 1 attempt .*Step 2: target 0.000000E[+]00.*Simulation completed successfully"
                   -DWORKING_DIR=${case_dir}
                   -DREMOVE_BEFORE=${case_dir}/simulation.sim
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_fepx_run_landing_np${num_ranks} PROPERTIES
                       TIMEOUT ${fepx_timeout_120}
                       PROCESSORS ${num_ranks}
                       LABELS mpi)
  set(case_dir ${CMAKE_CURRENT_BINARY_DIR}/app_runs/landing_max_bc_iter_np${num_ranks})
  file(MAKE_DIRECTORY ${case_dir})
  configure_file(app/cases/landing_max_bc_iter.cfg ${case_dir}/simulation.cfg COPYONLY)
  fepx_mesh_for_ranks(case_mesh regression/cases/meshes/sx2.msh ${num_ranks})
  configure_file(${case_mesh} ${case_dir}/simulation.msh COPYONLY)
  add_test(NAME test_fepx_run_landing_max_bc_iter_np${num_ranks}
           COMMAND ${CMAKE_COMMAND}
                   "-DCOMMAND=${command}"
                   -DEXPECTED_STATUS=nonzero
                   "-DEXPECTED_STDOUT=Error  : fepx: step 1, increment 3: no landing on the target  4.000000E[+]02 N within the tolerance  1.000E-01 in max_bc_iter 3 attempts [(]bracket: .*[)]. the last accepted increment is written. raise max_bc_iter or the tolerance"
                   -DEXPECTED_FILE=simulation.sim/results/elts/stress/stress.step1
                   -DWORKING_DIR=${case_dir}
                   -DREMOVE_BEFORE=${case_dir}/simulation.sim
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/core/check_run.cmake)
  set_tests_properties(test_fepx_run_landing_max_bc_iter_np${num_ranks} PROPERTIES
                       TIMEOUT ${fepx_timeout_120}
                       PROCESSORS ${num_ranks}
                       LABELS mpi)
endforeach()

# --- Development tools -----------------------------------------------------------------------

# Unit tests of tests/tools/compare_sim.py.
add_test(NAME tools_compare_sim
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/test_compare_sim.py)
set_tests_properties(tools_compare_sim PROPERTIES LABELS tools)

# Unit tests of dev/tools/check_quick_coverage.py (step 9.6; only where dev/ exists: it is on
# devel, not on main and not in the source tarballs).
if(EXISTS ${FEPX_ROOT_DIR}/dev/tools/check_quick_coverage.py)
  add_test(NAME tools_check_quick_coverage
           COMMAND ${Python3_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/tools/test_check_quick_coverage.py)
  set_tests_properties(tools_check_quick_coverage PROPERTIES LABELS tools)
endif()

# --- Verification and regression cases ------------------------------------------------------------

add_subdirectory(verification)
# The regression cases need the reference data (tests/regression/reference/), which the source
# tarball leaves out (step 9.8, ROADMAP D101: they are in the separate testdata tarball).
if(FEPX_HAVE_REFERENCES)
  add_subdirectory(regression)
else()
  message(STATUS "No regression reference data (tests/regression/reference/): the regression "
                 "cases are not defined (they come with the testdata tarball).")
endif()

# --- Manual (task M.1) -----------------------------------------------------------------------------

# The committed keyword pages (reference, compatibility table, every keyword's entry) must
# match the keyword table and the rules in the source (step 9.9, D104).
add_test(NAME docs_keywords_current
         COMMAND ${CMAKE_COMMAND}
                 -DGENERATOR=$<TARGET_FILE:fepx_keywords>
                 -DMANUAL_DIR=${FEPX_ROOT_DIR}/doc
                 -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/keyword_pages
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/docs/check_keyword_pages.cmake)
set_tests_properties(docs_keywords_current PROPERTIES LABELS docs)

# Every option-compatibility rule is used by the configuration builder and its message is
# expected by an error test (step 9.9, D104).
add_test(NAME docs_rules_used
         COMMAND ${CMAKE_COMMAND}
                 -DGENERATOR=$<TARGET_FILE:fepx_keywords>
                 -DPYTHON=${Python3_EXECUTABLE}
                 -DSCRIPT=${CMAKE_CURRENT_SOURCE_DIR}/docs/check_rules_used.py
                 -DSOURCE_DIR=${FEPX_ROOT_DIR}
                 -DRULES_FILE=${CMAKE_CURRENT_BINARY_DIR}/cfg_rules.txt
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/docs/check_rules_used.cmake)
set_tests_properties(docs_rules_used PROPERTIES LABELS docs)

# Every user-facing error and warning message has an entry in the manual's messages page
# (internal_error calls excepted; step 9.9, D104).
add_test(NAME docs_messages_documented
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/docs/check_messages.py
                 ${FEPX_ROOT_DIR})
set_tests_properties(docs_messages_documented PROPERTIES LABELS docs)

# Every configuration input shown in the manual is included from a tested file, not written in
# the page (step 9.9, D104).
add_test(NAME docs_manual_inputs
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/docs/check_manual_inputs.py
                 ${FEPX_ROOT_DIR}/doc)
set_tests_properties(docs_manual_inputs PROPERTIES LABELS docs)

# The committed slip-system lists must match crystal_mod's tables (step 6.1).
add_test(NAME docs_slip_systems_current
         COMMAND ${CMAKE_COMMAND}
                 -DGENERATOR=$<TARGET_FILE:fepx_slip_systems>
                 -DCOMMITTED=${FEPX_ROOT_DIR}/doc/slip_systems_table.rst
                 -DGENERATED=${CMAKE_CURRENT_BINARY_DIR}/slip_systems_table.rst
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/docs/check_keywords.cmake)
set_tests_properties(docs_slip_systems_current PROPERTIES LABELS docs)

# The tutorial pages' generated files exist (step 9.11): every file a page names, and the
# terminal output and tree of every case of tutorials/tutorials.txt
# (tutorials/tools/make_tutorial_pages.sh, run by hand, makes them).
add_test(NAME docs_tutorial_pages_current
         COMMAND ${Python3_EXECUTABLE} ${FEPX_ROOT_DIR}/tutorials/tools/tutorial_pages.py check
                 ${FEPX_ROOT_DIR}/tutorials/tutorials.txt
                 ${FEPX_ROOT_DIR}/doc)
set_tests_properties(docs_tutorial_pages_current PROPERTIES LABELS docs)

# The post-processing tutorial's script runs on the simple simulation's reference results (step
# 9.11; only where Matplotlib is installed and with the reference data).
execute_process(COMMAND ${Python3_EXECUTABLE} -c "import matplotlib"
                RESULT_VARIABLE fepx_matplotlib_missing OUTPUT_QUIET ERROR_QUIET)
if(NOT fepx_matplotlib_missing AND FEPX_HAVE_REFERENCES)
  add_test(NAME docs_tutorial_postprocess
           COMMAND ${Python3_EXECUTABLE} ${FEPX_ROOT_DIR}/tutorials/tools/tutorial_pages.py
                   postprocess-test
                   ${CMAKE_CURRENT_SOURCE_DIR}/regression/reference/fepx/manual_tutorial_simple_simulation
                   ${CMAKE_CURRENT_SOURCE_DIR}/regression/cases/meshes/tutorial_simple.msh
                   ${FEPX_ROOT_DIR}/tutorials/postprocessing/postprocess.py
                   ${CMAKE_CURRENT_BINARY_DIR}/tutorial_postprocess)
  set_tests_properties(docs_tutorial_postprocess PROPERTIES LABELS docs)
endif()

# The manual builds without warnings (only where Sphinx is installed; CI has a separate job). -E:
# a fresh environment every time (an incremental build reports a label moved to another page as
# a duplicate).
find_program(FEPX_SPHINX_BUILD sphinx-build)
if (FEPX_SPHINX_BUILD)
  add_test(NAME docs_manual_build
           COMMAND ${FEPX_SPHINX_BUILD} -E -W --keep-going -q -b html -d
                   ${CMAKE_CURRENT_BINARY_DIR}/manual_doctrees
                   ${FEPX_ROOT_DIR}/doc ${CMAKE_CURRENT_BINARY_DIR}/manual_html)
  set_tests_properties(docs_manual_build PROPERTIES LABELS docs)
endif()

# The test programs that distribute a mesh on any number of ranks (D108).
foreach(target test_sim_writer test_sim_writer_errors test_mesh_local test_mesh_distribute test_halo
               test_constraint test_fe_errors test_rigid_rotation test_pcg_precond test_periodic
               test_increment test_fe_solve test_lame test_isotropic_start)
  target_link_libraries(${target} PRIVATE fepx_test_mesh)
endforeach()
