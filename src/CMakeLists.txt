# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# FEPX top-level build.
#
# Usage (see README.md and CLAUDE.md §7):
#   cmake -S src -B build -DCMAKE_BUILD_TYPE=Release        # or Debug, or Checked
#   cmake --build build -j
#   ctest --test-dir build --output-on-failure
#
# Options:
#   FEPX_WERROR   treat compiler warnings as errors (ON in CI)
#   FEPX_INSTALL_RPATH  the installed executable keeps the run path of its libraries (default ON;
#                  OFF for a package, e.g. conda's, step 9.8)
#   BUILD_TESTING  build the test suite (default ON)
#   FEPX_ENABLE_HDF5  the optional HDF5 output (output_format hdf5, step 9.1, ROADMAP D93):
#                  AUTO (the default, step 9.8, D101: built when found), ON (required) or OFF;
#                  needs the serial HDF5 library with its Fortran interface (e.g. Debian/Ubuntu
#                  libhdf5-dev), and h5ls, h5dump, h5diff (hdf5-tools) for the tests

cmake_minimum_required(VERSION 3.20)

# Checked is not one of CMake's own build types: the initial value of its optimization flags
# (the run-time checks are added in src/cmake/CompilerFlags.cmake) must be given before project(),
# which creates the cache entry CMAKE_Fortran_FLAGS_CHECKED from it (set after project(), the
# entry already exists, empty, and the build is not optimized).
set(CMAKE_Fortran_FLAGS_CHECKED_INIT "-O2 -g")

project(fepx VERSION 3.0.0 LANGUAGES Fortran)

# The repository's root: the build's source directory is src/ (FEPX-dev's layout, step 9.21,
# ROADMAP D118); tests/, doc/ and tutorials/ are next to it.
get_filename_component(FEPX_ROOT_DIR ${PROJECT_SOURCE_DIR}/.. ABSOLUTE)

# Version suffix (step 9.7, Q58): "-dev" until the release (step 9.10) sets the release version.
set(FEPX_VERSION_SUFFIX "")

# The version and the git tag (step 9.8, ROADMAP D101): the version above is the only source of
# the version (a source tarball has no git); a commit tagged vX.Y.Z must have the version X.Y.Z
# and no suffix. Checked only in the source tree's own git repository.
find_package(Git QUIET)
if(GIT_FOUND)
  execute_process(COMMAND ${GIT_EXECUTABLE} -C ${FEPX_ROOT_DIR} rev-parse --show-toplevel
                  OUTPUT_VARIABLE fepx_git_top RESULT_VARIABLE fepx_git_status
                  OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(fepx_git_status EQUAL 0)
    get_filename_component(fepx_git_top "${fepx_git_top}" REALPATH)
    get_filename_component(fepx_source_real "${FEPX_ROOT_DIR}" REALPATH)
    if(fepx_git_top STREQUAL fepx_source_real)
      # The version tags of HEAD (none on an untagged commit).
      execute_process(COMMAND ${GIT_EXECUTABLE} -C ${FEPX_ROOT_DIR} tag --points-at HEAD
                              --list "v[0-9]*"
                      OUTPUT_VARIABLE fepx_version_tags OUTPUT_STRIP_TRAILING_WHITESPACE
                      ERROR_QUIET)
      string(REPLACE "\n" ";" fepx_version_tags "${fepx_version_tags}")
      foreach(tag IN LISTS fepx_version_tags)
        if(NOT tag STREQUAL "v${PROJECT_VERSION}${FEPX_VERSION_SUFFIX}"
           OR NOT FEPX_VERSION_SUFFIX STREQUAL "")
          message(FATAL_ERROR "The commit is tagged ${tag}, but the version in CMakeLists.txt is "
                              "${PROJECT_VERSION}${FEPX_VERSION_SUFFIX}: a release tag vX.Y.Z "
                              "needs project(VERSION X.Y.Z) and an empty FEPX_VERSION_SUFFIX.")
        endif()
      endforeach()
    endif()
  endif()
endif()

# Build types (single-config generators):
#   Release  optimized (CMake's default -O3), no run-time checks; the default
#   Debug    -O0 -g and run-time checks (cmake/CompilerFlags.cmake), for interactive debugging
#   Checked  -O2 -g and the same run-time checks as Debug; the CI's checked build (ROADMAP D88:
#            measured 1.2-1.6 times Release on the simulation tests, Debug 4-6 times)
# Default to a Release build when no build type is given (single-config generators only).
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type: Debug, Checked or Release" FORCE)
endif()
if(NOT CMAKE_CONFIGURATION_TYPES)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Checked Release)
  if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Checked|Release)$")
    message(FATAL_ERROR "CMAKE_BUILD_TYPE '${CMAKE_BUILD_TYPE}' is not one of Debug, Checked, Release.")
  endif()
  # A Checked build without optimization would be a slow Debug build (a build directory
  # configured before CMAKE_Fortran_FLAGS_CHECKED_INIT existed has the entry empty).
  if(CMAKE_BUILD_TYPE STREQUAL "Checked" AND NOT CMAKE_Fortran_FLAGS_CHECKED MATCHES "-O[1-3s]")
    message(FATAL_ERROR "CMAKE_Fortran_FLAGS_CHECKED is '${CMAKE_Fortran_FLAGS_CHECKED}', without "
                        "optimization: configure a new build directory, or give "
                        "-DCMAKE_Fortran_FLAGS_CHECKED=\"-O2 -g\".")
  endif()
endif()

option(FEPX_WERROR "Treat compiler warnings as errors" OFF)
option(BUILD_TESTING "Build the FEPX test suite" ON)
set(FEPX_ENABLE_HDF5 AUTO CACHE STRING "Build the optional HDF5 output (output_format hdf5): AUTO, ON or OFF")
set_property(CACHE FEPX_ENABLE_HDF5 PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${FEPX_ENABLE_HDF5}" fepx_hdf5_request)
if(NOT fepx_hdf5_request MATCHES "^(AUTO|ON|OFF)$")
  message(FATAL_ERROR "FEPX_ENABLE_HDF5 '${FEPX_ENABLE_HDF5}' is not one of AUTO, ON, OFF.")
endif()

# MPI: any implementation is accepted, but it must provide the Fortran 2008 bindings
# (`use mpi_f08`), which are the only MPI interface FEPX uses (ROADMAP D5).
# MPI_DETERMINE_LIBRARY_VERSION makes FindMPI report the implementation (Open MPI, MPICH, ...)
# in MPI_Fortran_LIBRARY_VERSION_STRING; the test helpers use it to choose launcher flags.
set(MPI_DETERMINE_LIBRARY_VERSION TRUE)
find_package(MPI REQUIRED COMPONENTS Fortran)
if(NOT MPI_Fortran_HAVE_F08_MODULE)
  message(FATAL_ERROR
    "The MPI installation found (${MPI_Fortran_LIBRARY_VERSION_STRING}) does not provide "
    "the mpi_f08 Fortran module, which FEPX requires.")
endif()

# HDF5 (optional, step 9.1, ROADMAP D10, D93; AUTO since step 9.8, D101): the serial library,
# written by rank 0 only, so any MPI may be used with it. A parallel HDF5 is accepted with a
# warning: it works when it was built with the same MPI as FEPX (it is used serially).
# FEPX_HAVE_HDF5 (ON or OFF) is the result, used by the rest of the build: with ON the library
# (and, with the test suite, the HDF5 tools) must be found; with AUTO the output is built when
# they are found.
set(FEPX_HAVE_HDF5 OFF)
if(NOT fepx_hdf5_request STREQUAL "OFF")
  set(HDF5_PREFER_PARALLEL FALSE)
  # An HDF5 built with CMake (always the case from HDF5 2.0) installs a CMake package file
  # (hdf5-config.cmake), which FindHDF5 uses first; for a thread-safe HDF5 (e.g. conda-forge's),
  # that file looks for the threads library (FindThreads), which works only with the C language
  # enabled, and FEPX's project is Fortran only. C is therefore enabled when a C compiler is
  # available (it comes with gfortran, gcc, or with ifx, icx); otherwise the package file is not
  # used and FindHDF5 looks for HDF5 through its compiler wrapper (h5fc), as for an HDF5 built
  # without CMake (e.g. Debian's libhdf5-dev, HDF5 1.x).
  include(CheckLanguage)
  check_language(C)
  if(CMAKE_C_COMPILER)
    enable_language(C)
  else()
    set(HDF5_NO_FIND_PACKAGE_CONFIG_FILE ON)
  endif()
  if(fepx_hdf5_request STREQUAL "ON")
    find_package(HDF5 REQUIRED COMPONENTS Fortran)
  else()
    find_package(HDF5 QUIET COMPONENTS Fortran)
  endif()
  if(HDF5_FOUND)
    set(FEPX_HAVE_HDF5 ON)
    if(HDF5_IS_PARALLEL)
      message(WARNING "The HDF5 found (${HDF5_Fortran_LIBRARIES}) is a parallel build; FEPX "
                      "uses it serially, which needs the same MPI as FEPX. A serial HDF5 "
                      "(e.g. libhdf5-dev) is recommended.")
    endif()
  endif()
  # The HDF5 tests read the files with the HDF5 tools (hdf5-tools: h5ls, h5dump, h5diff).
  if(FEPX_HAVE_HDF5 AND BUILD_TESTING)
    set(missing_tools "")
    foreach(tool h5ls h5dump h5diff)
      find_program(FEPX_${tool}_PROGRAM ${tool})
      if(NOT FEPX_${tool}_PROGRAM)
        list(APPEND missing_tools ${tool})
      endif()
    endforeach()
    if(missing_tools AND fepx_hdf5_request STREQUAL "ON")
      message(FATAL_ERROR "The HDF5 tests need ${missing_tools} (e.g. apt install hdf5-tools), "
                          "or configure with -DBUILD_TESTING=OFF.")
    elseif(missing_tools)
      message(WARNING "HDF5 was found but not ${missing_tools}, which the HDF5 tests need "
                      "(e.g. apt install hdf5-tools): the HDF5 output is not built. Install "
                      "them, or configure with -DFEPX_ENABLE_HDF5=ON -DBUILD_TESTING=OFF.")
      set(FEPX_HAVE_HDF5 OFF)
    endif()
    set(FEPX_H5DIFF ${FEPX_h5diff_PROGRAM})
  endif()
endif()

# All Fortran module files go to one directory in the build tree.
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)

# Compiler flag sets (gfortran and ifx).
include(cmake/CompilerFlags.cmake)

# Installation layout (step 9.8, ROADMAP D101): GNUInstallDirs under CMAKE_INSTALL_PREFIX; the
# package (and the documentation directory share/doc/<package>) is named fepx (the project keeps
# its name fepx until the renaming of step 9.10). With FEPX_INSTALL_RPATH (the default), the
# installed executable keeps the run path of the libraries it was linked with outside the
# system directories (an MPI or an HDF5 of a module or /opt), so that it runs without
# LD_LIBRARY_PATH; a package is built without it (e.g. conda sets the run paths of its packages
# itself).
set(FEPX_PACKAGE_NAME fepx)
if(NOT DEFINED CMAKE_INSTALL_DOCDIR)
  set(CMAKE_INSTALL_DOCDIR share/doc/${FEPX_PACKAGE_NAME} CACHE PATH
      "documentation root (DATAROOTDIR/doc/${FEPX_PACKAGE_NAME})")
endif()
include(GNUInstallDirs)
option(FEPX_INSTALL_RPATH "Keep the run path of the linked libraries in the installed executable" ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ${FEPX_INSTALL_RPATH})

# FEPX sources (library fepx_core and the executables; the former src/CMakeLists.txt, merged
# here when the build's source directory became src/, step 9.21, D118).
#
# Layout of the sources (this directory):
#   core/  kinds, mathematical constants, fatal-error handling (step 1.1), logging (2.1), the
#          build information (9.7, generated from build_info_mod.f90.in)
#   math/  Mandel tensor algebra (step 1.2), frame transforms and orientation descriptors (1.3),
#          small dense linear systems (4.1), time integration of a state vector (6.4)
#   material/  crystal structures, slip systems, crystal-frame stiffness (1.4); material state,
#              the material-solver contract and the elastic-only model (4.1); slip families
#              and the slip kinetics, power law and [M*], [M*_tan] (6.1)
#   element/  element library, one module per element with its quadrature rule (D36):
#             10-node tetrahedron (step 3.1), 20-node hexahedron (3.2), face elements (3.3);
#             element_mod: type-independent descriptor and geometry kernels (3.4)
#   fe/    finite element formulation: element matrices, the step schedule (5.1) and the
#          pure parts of the boundary-velocity iteration (8.2b) (in fepx_core), element-by-element
#          operations, Jacobi-PCG and Dirichlet conditions (in fepx_parallel) (step 4.2); face-set
#          forces and the increment solve (successive approximation) (in fepx_parallel) (4.4);
#          velocity constraints v = T v_hat + g, multi-point constraints (8.3a, in fepx_parallel)
#   io/    text input files, Neper mesh file parser (2.2); configuration keyword table, reader
#          and typed configuration (4.3); .sim output writer (4.3, in fepx_parallel); HDF5
#          files (9.1: hdf5_file_mod.f90 with FEPX_ENABLE_HDF5, otherwise its stub
#          hdf5_file_stub.f90, the same module whose routines are fatal errors)
#   mesh/  rank-local mesh type (2.2); its distribution to the ranks (2.2), local numbering and
#          communication graph (2.3) (the last two in fepx_parallel)
#   app/   the fepx executable (step 4.4)
#   parallel/  MPI environment and collectives (2.1), halo exchange and distributed vector
#              operations (2.4), the transfer of the constraint columns (8.3a, D74): the only
#              code that uses MPI (conventions §6)
#
# The MPI-free modules are compiled into the static library fepx_core, which every test
# executable links; src/parallel/ is the library fepx_parallel (fepx_core + MPI), linked by the
# MPI tests (and the fepx executable).

# The HDF5 file module, or its stub (step 9.1).

if(FEPX_HAVE_HDF5)
  set(fepx_hdf5_source io/hdf5_file_mod.f90)
else()
  set(fepx_hdf5_source io/hdf5_file_stub.f90)
endif()

# The build-information module (step 9.7): generated from its template at configure time and
# again at every build (version, git hash and commit date, dirty flag), rewritten only when its
# text changes (cmake/BuildInfo.cmake).
set(fepx_build_info_source ${CMAKE_CURRENT_BINARY_DIR}/build_info_mod.f90)
set(fepx_build_info_command ${CMAKE_COMMAND} -DSOURCE_DIR=${FEPX_ROOT_DIR}
    -DTEMPLATE=${CMAKE_CURRENT_SOURCE_DIR}/core/build_info_mod.f90.in
    -DOUTPUT=${fepx_build_info_source} -DVERSION=${PROJECT_VERSION} -DSUFFIX=${FEPX_VERSION_SUFFIX}
    -P ${PROJECT_SOURCE_DIR}/cmake/BuildInfo.cmake)
execute_process(COMMAND ${fepx_build_info_command})
add_custom_target(fepx_build_info
                  COMMAND ${fepx_build_info_command}
                  BYPRODUCTS ${fepx_build_info_source}
                  COMMENT "Build information (version, git hash and date)")

add_library(fepx_core STATIC
  core/kinds_mod.f90
  ${fepx_build_info_source}
  core/constants_mod.f90
  core/error_mod.f90
  core/log_mod.f90
  math/tensor_mod.f90
  math/rotation_mod.f90
  math/orientation_mod.f90
  math/linalg_mod.f90
  math/time_integration_mod.f90
  material/crystal_mod.f90
  material/slip_kinetics_mod.f90
  material/hardening_mod.f90
  material/material_state_mod.f90
  material/material_model_mod.f90
  material/elastic_model_mod.f90
  material/yield_vertices_mod.f90
  material/crystal_plasticity_model_mod.f90
  element/tet10_mod.f90
  element/hex20_mod.f90
  element/element_mod.f90
  fe/element_matrix_mod.f90
  io/text_file_mod.f90
  io/cfg_table_mod.f90
  io/cfg_rules_mod.f90
  io/cfg_manual_mod.f90
  io/cfg_store_mod.f90
  io/config_mod.f90
  fe/step_schedule_mod.f90
  fe/bc_iteration_mod.f90
  io/restart_file_mod.f90
  io/mesh_file_mod.f90
  ${fepx_hdf5_source}
  mesh/mesh_mod.f90
)

add_dependencies(fepx_core fepx_build_info)

# HDF5 (step 9.1): the library and its Fortran modules; the stub's arguments are unused.
if(FEPX_HAVE_HDF5)
  target_include_directories(fepx_core PUBLIC ${HDF5_Fortran_INCLUDE_DIRS})
  target_link_libraries(fepx_core PUBLIC ${HDF5_Fortran_LIBRARIES})
else()
  set_source_files_properties(io/hdf5_file_stub.f90 PROPERTIES
                              COMPILE_OPTIONS "$<$<Fortran_COMPILER_ID:GNU>:-Wno-unused-dummy-argument>")
endif()

add_library(fepx_parallel STATIC
  parallel/parallel_mod.f90
  parallel/halo_mod.f90
  mesh/mesh_local_mod.f90
  mesh/mesh_distribute_mod.f90
  fe/ebe_mod.f90
  parallel/constraint_exchange_mod.f90
  fe/constraint_mod.f90
  fe/periodic_mod.f90
  fe/pcg_mod.f90
  fe/bc_mod.f90
  fe/face_force_mod.f90
  fe/increment_mod.f90
  fe/isotropic_start_mod.f90
  io/sim_writer_mod.f90
)
target_link_libraries(fepx_parallel PUBLIC fepx_core MPI::MPI_Fortran)

# The FEPX executable (step 4.4: one elastic increment).
add_executable(fepx app/fepx.f90)
target_link_libraries(fepx PRIVATE fepx_parallel)

# The keyword reference generator of the manual (task M.1; serial).
add_executable(fepx_keywords app/fepx_keywords.f90)
target_link_libraries(fepx_keywords PRIVATE fepx_core)

# The slip-system lists of the manual (step 6.1; serial).
add_executable(fepx_slip_systems app/fepx_slip_systems.f90)
target_link_libraries(fepx_slip_systems PRIVATE fepx_core)

# Installation (step 9.8, ROADMAP D101): the fepx executable only; fepx_keywords and
# fepx_slip_systems make the manual's tables and are not installed, nor are the libraries.
install(TARGETS fepx RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Installation and the source archives (step 9.8).
include(cmake/Packaging.cmake)

# The conda packages (recipe in conda/recipe, one package per MPI implementation), built from the
# local source tree and tested with rattler-build, in conda/<platform>/ of the build directory.
find_program(RATTLER_BUILD rattler-build)
mark_as_advanced(RATTLER_BUILD)
if(RATTLER_BUILD)
  add_custom_target(conda
    COMMAND ${RATTLER_BUILD} build -r ${FEPX_ROOT_DIR}/conda/recipe
            -c conda-forge --output-dir ${PROJECT_BINARY_DIR}/conda
    WORKING_DIRECTORY ${FEPX_ROOT_DIR}
    COMMENT "Building and testing the conda packages in ${PROJECT_BINARY_DIR}/conda"
    USES_TERMINAL
    VERBATIM)
else()
  add_custom_target(conda
    COMMAND ${CMAKE_COMMAND} -E echo "rattler-build not found: install it (e.g., conda install"
            "-n base -c conda-forge rattler-build) and rerun cmake."
    COMMAND ${CMAKE_COMMAND} -E false
    VERBATIM)
endif()

# Test suite.
if(BUILD_TESTING)
  # Python 3 + NumPy: development tool for comparing outputs in tests (CLAUDE.md §6).
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy"
                  RESULT_VARIABLE numpy_status OUTPUT_QUIET ERROR_QUIET)
  if(NOT numpy_status EQUAL 0)
    message(FATAL_ERROR "The test suite needs NumPy for ${Python3_EXECUTABLE} "
                        "(e.g. apt install python3-numpy), or configure with -DBUILD_TESTING=OFF.")
  endif()
  # The regression reference data (absent from the source tarball, step 9.8, D101).
  if(IS_DIRECTORY ${FEPX_ROOT_DIR}/tests/regression/reference/fepx)
    set(FEPX_HAVE_REFERENCES ON)
  else()
    set(FEPX_HAVE_REFERENCES OFF)
  endif()
  enable_testing()
  include(cmake/FEPXTesting.cmake)
  add_subdirectory(${FEPX_ROOT_DIR}/tests ${PROJECT_BINARY_DIR}/tests)
endif()

message(STATUS "FEPX ${PROJECT_VERSION}: build type ${CMAKE_BUILD_TYPE}, "
               "compiler ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}, "
               "MPI: ${MPI_Fortran_LIBRARY_VERSION_STRING}, "
               "HDF5 output: ${FEPX_HAVE_HDF5}")
