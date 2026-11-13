# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# CTest helpers for FEPX.
#
#   fepx_add_unit_test(<name> SOURCES <file>...)
#       Serial test executable <name>, linked with the assertion module and the FEPX library
#       fepx_core; one CTest test <name>.
#
#   fepx_add_mpi_test(<name> SOURCES <file>... RANKS <n>... [ARGS <arg>...])
#       MPI test executable <name>, linked with the assertion module and fepx_parallel (which
#       brings fepx_core and MPI); one CTest test <name>_np<n> per rank
#       count, run through MPIEXEC_EXECUTABLE with the given command-line ARGS. CLAUDE.md
#       requires 1, 2 and 4 ranks for MPI tests.
#
#   fepx_mpiexec_command(<out_var> <num_ranks> <executable>)
#       The launcher command (a CMake list) that runs <executable> on <num_ranks> ranks, as used
#       by fepx_add_mpi_test; for tests that wrap the run in a script (e.g. check_run.cmake).
#
#   fepx_read_quick_tests(<file>)
#       Reads the names of the quick tests other than the unit tests (tests/quick_tests.txt: one
#       test name per line, # comments) into the global property FEPX_QUICK_TESTS.
#
#   fepx_defer_test_tiers()
#       Called at the top of every CMakeLists.txt that defines tests: schedules
#       fepx_apply_test_tiers for the end of the directory (cmake_language(DEFER)), which runs
#       even when the file ends early with return() (a call at the end of the file was skipped
#       there and left its tests without tier labels, found 2026-10-24).
#
#   fepx_apply_test_tiers()
#       Adds the tier labels (ROADMAP D88) to the tests of the current directory: "full" to
#       every test but the tutorial cases (label "tutorial"), "quick" to every unit test and to
#       every test named in FEPX_QUICK_TESTS. The names labelled are recorded in
#       FEPX_QUICK_TESTS_FOUND, the directory in FEPX_TIERED_DIRS.
#
#   fepx_check_quick_tests()
#       Deferred at the top of tests/CMakeLists.txt after its own fepx_defer_test_tiers, so
#       that it runs after every test directory. Configuration errors: a directory under tests/
#       with tests whose tiers were not applied (every test must be full or tutorial); a name in
#       FEPX_QUICK_TESTS that defines no test (a renamed or removed test must leave the list),
#       except docs_manual_build without Sphinx; a quick tutorial case.
#
# Every test gets a label: "unit" or "mpi" (select with `ctest -L mpi`, exclude with `-LE mpi`),
# and the tier labels "quick", "full" (fepx_apply_test_tiers; select with `ctest -L quick`).

# Open MPI launcher flags (none are needed for other implementations):
#   --oversubscribe: Open MPI refuses to start more ranks than available cores unless asked to;
#     CI runners and laptops may have fewer cores than the largest rank count tested.
#   --mca btl_tcp_if_include lo: the tests run on one machine, so Open MPI's TCP transport is
#     restricted to the loopback interface. Otherwise it tries to pair every network interface,
#     which on the development laptop (Wi-Fi with changing IPv6 addresses, an address-less
#     interface) occasionally failed at startup with "Unable to find reachable pairing between
#     local and remote interfaces" under `ctest -j 8` (ROADMAP Q11). The loopback interface is
#     `lo` on Linux and `lo0` on macOS; with `lo` on macOS, Open MPI ignores the value and warns
#     on every run, and the warning can put a blank line in the run's standard output
#     (ROADMAP D112).
# FEPX_MPIEXEC_EXTRA_FLAGS is a CMake list (for commands built here);
# FEPX_MPIEXEC_EXTRA_FLAGS_TEXT is the same, space-separated, for the -DMPIEXEC_EXTRA=...
# arguments of the run scripts (split again with separate_arguments).
set(FEPX_MPIEXEC_EXTRA_FLAGS "")
if(MPI_Fortran_LIBRARY_VERSION_STRING MATCHES "Open MPI")
  # The loopback interface's name.
  set(loopback lo)
  if(APPLE)
    set(loopback lo0)
  endif()
  set(FEPX_MPIEXEC_EXTRA_FLAGS --oversubscribe --mca btl_tcp_if_include ${loopback})
endif()
string(JOIN " " FEPX_MPIEXEC_EXTRA_FLAGS_TEXT ${FEPX_MPIEXEC_EXTRA_FLAGS})

function(fepx_add_unit_test name)
  cmake_parse_arguments(ARG "" "" "SOURCES" ${ARGN})
  add_executable(${name} ${ARG_SOURCES})
  target_link_libraries(${name} PRIVATE fepx_testing fepx_core)
  add_test(NAME ${name} COMMAND ${name})
  set_tests_properties(${name} PROPERTIES LABELS unit)
endfunction()

function(fepx_mpiexec_command out_var num_ranks executable)
  set(command ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${num_ranks}
              ${FEPX_MPIEXEC_EXTRA_FLAGS} ${MPIEXEC_PREFLAGS} ${executable} ${MPIEXEC_POSTFLAGS})
  set(${out_var} "${command}" PARENT_SCOPE)
endfunction()

# The file of a test mesh for a number of ranks (ROADMAP D108): FEPX runs a mesh on as many ranks
# as it has partitions, so a mesh that runs on another number of ranks has a copy with that many
# partitions: <mesh>_part<n>.msh (Neper's table). Returns the copy when it exists, else the mesh
# itself; 1 rank takes any mesh. A relative path is relative to the current source directory.
function(fepx_mesh_for_ranks out mesh num_ranks)
  set(${out} "${mesh}" PARENT_SCOPE)
  if(num_ranks EQUAL 1)
    return()
  endif()
  string(REGEX REPLACE "(_part[0-9]+)?[.]msh$" "" stem "${mesh}")
  set(candidate "${stem}_part${num_ranks}.msh")
  get_filename_component(candidate_abs "${candidate}" ABSOLUTE BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
  if(EXISTS "${candidate_abs}")
    set(${out} "${candidate}" PARENT_SCOPE)
  endif()
endfunction()

function(fepx_add_mpi_test name)
  cmake_parse_arguments(ARG "" "" "SOURCES;RANKS;ARGS" ${ARGN})
  add_executable(${name} ${ARG_SOURCES})
  target_link_libraries(${name} PRIVATE fepx_testing fepx_parallel)
  foreach(num_ranks IN LISTS ARG_RANKS)
    fepx_mpiexec_command(command ${num_ranks} $<TARGET_FILE:${name}>)
    add_test(NAME ${name}_np${num_ranks} COMMAND ${command} ${ARG_ARGS})
    # FEPX_EXPECTED_RANKS lets test_parallel check that the launcher started one job of
    # num_ranks ranks (not num_ranks singletons).
    set_tests_properties(${name}_np${num_ranks} PROPERTIES
                         PROCESSORS ${num_ranks}
                         ENVIRONMENT FEPX_EXPECTED_RANKS=${num_ranks}
                         LABELS mpi)
  endforeach()
endfunction()

function(fepx_read_quick_tests file)
  # Each line without its comment and surrounding blanks; empty lines skipped.
  file(STRINGS ${file} lines)
  set(names "")
  # (A line "hdf5: <name>" names a test of the HDF5 output, defined only with FEPX_ENABLE_HDF5,
  # step 9.1; it is skipped in a build without it. A line "dev: <name>" names a test of a tool
  # of dev/, defined only where dev/ exists: it is skipped on main and in the source tarballs.)
  foreach(line IN LISTS lines)
    string(REGEX REPLACE "#.*$" "" line "${line}")
    string(STRIP "${line}" line)
    if(line MATCHES "^hdf5:[ ]*(.+)$")
      if(NOT FEPX_HAVE_HDF5)
        continue()
      endif()
      set(line ${CMAKE_MATCH_1})
    endif()
    if(line MATCHES "^dev:[ ]*(.+)$")
      if(NOT EXISTS ${FEPX_ROOT_DIR}/dev)
        continue()
      endif()
      set(line ${CMAKE_MATCH_1})
    endif()
    if(NOT line STREQUAL "")
      list(APPEND names ${line})
    endif()
  endforeach()
  set_property(GLOBAL PROPERTY FEPX_QUICK_TESTS ${names})
  set_property(GLOBAL PROPERTY FEPX_QUICK_TESTS_FOUND "")
  # Re-run the configuration when the list changes.
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${file})
endfunction()

function(fepx_defer_test_tiers)
  cmake_language(DEFER CALL fepx_apply_test_tiers)
endfunction()

function(fepx_apply_test_tiers)
  set_property(GLOBAL APPEND PROPERTY FEPX_TIERED_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
  get_property(quick_names GLOBAL PROPERTY FEPX_QUICK_TESTS)
  get_property(tests DIRECTORY PROPERTY TESTS)
  foreach(test IN LISTS tests)
    get_test_property(${test} LABELS labels)
    if(NOT labels)
      set(labels "")
    endif()
    # Every test but the tutorial cases is in the full tier.
    if(NOT "tutorial" IN_LIST labels)
      list(APPEND labels full)
    endif()
    # The unit tests and the listed tests are in the quick tier.
    if("unit" IN_LIST labels OR test IN_LIST quick_names)
      list(APPEND labels quick)
      set_property(GLOBAL APPEND PROPERTY FEPX_QUICK_TESTS_FOUND ${test})
      if("tutorial" IN_LIST labels)
        message(FATAL_ERROR "tests/quick_tests.txt: ${test} is a tutorial case, not a quick test.")
      endif()
    endif()
    set_tests_properties(${test} PROPERTIES LABELS "${labels}")
  endforeach()
endfunction()

function(fepx_check_quick_tests)
  get_property(quick_names GLOBAL PROPERTY FEPX_QUICK_TESTS)
  get_property(found GLOBAL PROPERTY FEPX_QUICK_TESTS_FOUND)
  get_property(tiered_dirs GLOBAL PROPERTY FEPX_TIERED_DIRS)
  # Every directory under this one (tests/) that defines tests must have had its tiers applied.
  set(pending_dirs ${CMAKE_CURRENT_SOURCE_DIR})
  while(pending_dirs)
    list(POP_FRONT pending_dirs dir)
    get_property(dir_tests DIRECTORY ${dir} PROPERTY TESTS)
    if(dir_tests AND NOT dir IN_LIST tiered_dirs)
      message(FATAL_ERROR "${dir}/CMakeLists.txt defines tests without tier labels: call "
                          "fepx_defer_test_tiers() at its top (cmake/FEPXTesting.cmake).")
    endif()
    get_property(subdirs DIRECTORY ${dir} PROPERTY SUBDIRECTORIES)
    list(APPEND pending_dirs ${subdirs})
  endwhile()
  # Every listed name must define a test.
  foreach(name IN LISTS quick_names)
    if(name IN_LIST found)
      continue()
    endif()
    # Tests that are defined only in some configurations.
    if(NOT FEPX_SPHINX_BUILD AND name STREQUAL "docs_manual_build")
      continue()
    endif()
    if(NOT FEPX_HAVE_REFERENCES AND name MATCHES "^regression_")
      continue()
    endif()
    message(FATAL_ERROR "tests/quick_tests.txt: no test named ${name}.")
  endforeach()
endfunction()
