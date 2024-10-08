# This file is part of the FEPX software package.
# Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Top-level CMakeLists.txt for building FEPX
cmake_minimum_required(VERSION 3.0)
project(fepx Fortran)
set(FEPX_VERSION \"2.1.0\")

# Prepare CMake for Fortran and locate MPI package
enable_language(Fortran)
find_package(MPI REQUIRED)

# Set the build type for appropriate compile flags
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build: Debug Release" FORCE)
endif()

set(BUILD_TESTING_MODE "Normal" CACHE STRING "Choose the testing mode: Normal Minimal Writing")
set(BUILD_TESTING_DIFF "Soft" CACHE STRING "Choose the testing method: Soft Hard")

set(OUTPUT_PRECISION "1e-12" CACHE STRING "Output precision (for real numbers)")

# Option for verbose build
option(VERBOSE "Build verbose Makefiles" OFF)
if(VERBOSE)
    set(CMAKE_VERBOSE_MAKEFILE ON)
else(VERBOSE)
    set(CMAKE_VERBOSE_MAKEFILE OFF)
endif(VERBOSE)

# Option to suppress IEEE warnings on program exit
option(IEEE_SUPPRESS "Suppress IEEE signals on program exit" ON)
if(IEEE_SUPPRESS)
    set(CMAKE_Fortran_FLAGS "-ffpe-summary='none'")
else(IEEE_SUPPRESS)
    set(CMAKE_Fortran_FLAGS "")
endif(IEEE_SUPPRESS)

# Code to silence warnings from RANLIB on OSX
# From https://stackoverflow.com/questions/4929255/building-static-libraries-on-mac-using-cmake-and-gcc
if (APPLE)
  SET(CMAKE_Fortran_ARCHIVE_CREATE   "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
  SET(CMAKE_Fortran_ARCHIVE_FINISH   "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
endif ()

# Check that this machine supports F90 compilation
if(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    message(FATAL_ERROR "Fortran compiler does not support F90")
endif(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

if((${CMAKE_Fortran_COMPILER_ID} MATCHES Intel) OR (${CMAKE_Fortran_COMPILER_ID} MATCHES GNU))
    message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}, ${CMAKE_Fortran_COMPILER_ID}")
else()
    message(FATAL_ERROR "Unsupported Fortran compiler (use Intel or GNU).")
endif()

# Set compiler flag options depending on the compiler located
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_MODDIR_FLAG   "-J")
#    set(CMAKE_Fortran_FLAGS         "")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fcheck=all -fbacktrace -fbounds-check -Wall -Wno-argument-temporarily-ignored -Wno-unused-dummy-argument")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ftree-vectorize -funroll-loops")
# IFORT is currently not supported, but will be in the future
#ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
#   set(CMAKE_Fortran_MODDIR_FLAG   "")
#   set(CMAKE_Fortran_FLAGS         "")
#   set(CMAKE_Fortran_FLAGS_DEBUG   "")
#   set(CMAKE_Fortran_FLAGS_RELEASE "")
endif()

# Set directories for CMake to access files
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# Preamble info
STRING(TIMESTAMP CONFIG_DATE "%Y-%m-%d")
message(STATUS "Configuring Makefiles for the FEPX software package...")
message(STATUS "Built date: (${CONFIG_DATE})")
message(STATUS "OS: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}") # Need to grab input cache here?
message(STATUS "Version: ${VERSION}")

# configure a header file to pass some of the CMake settings to the
# source code
configure_file(
  "${PROJECT_SOURCE_DIR}/fepx_config.f90.in"
  "${PROJECT_BINARY_DIR}/fepx_config.f90"
)

# Add in the subdirectories and includes

file(GLOB fepx
    ${PROJECT_BINARY_DIR}/fepx_config.f90
    fepx.f90
    fepx2.f90
    
    read_input/read_input.f90
    read_input/read_input_cfg.f90
    read_input/read_input_msh.f90
    read_input/read_input_msh2.f90
    read_input/boundary_conditions_mod1.f90
    read_input/boundary_conditions_mod2.f90
    read_input/boundary_conditions_mod3.f90
    read_input/boundary_conditions_mod2_legacy.f90
    read_input/read_input_mpcs/multi_point_constraints_mod1.f90
    read_input/read_input_mpcs/multi_point_constraints_mod2.f90
    
    read_input/read_input_periodicity/periodicity_mod1.f90
    read_input/read_input_periodicity/periodicity_mod2.f90
    read_input/read_input_periodicity/periodicity_mod3.f90
    
    general/crys_type_mod1.f90
    general/crys_type_mod2.f90
    general/crys_type_mod3.f90
    general/general_mod.f90
    general/types_mod.f90
    general/printing_type_mod.f90
    general/exec_type_mod.f90
    general/loading_type_mod1.f90
    general/loading_type_mod2.f90
    general/loading_options_type_mod.f90
    general/mesh_type_mod.f90
    general/utils.f90
    general/orientation_conversion_mod.f90
    parallel/gather_scatter_mod.f90
    parallel/parallel_mod.f90
    res_init/res_init_mod1.f90
    res_init/res_init_mod2.f90

    triaxial/driver_triaxclr_mod.f90
    triaxial/driver_triaxcsr_mod.f90
    triaxial/driver_triax_utilities_mod.f90

    libfepx/driver_uniaxial_control_mod1.f90
    libfepx/driver_uniaxial_control_mod2.f90

    libfepx/solveit_isovp/solveit_isovp_mod1.f90
    libfepx/solveit_isovp/solveit_isovp_mod2.f90
    libfepx/solveit_isovp/solveit_isovp_mod3.f90
    libfepx/solveit_isovp/solveit_isovp_mod4.f90

    libfepx/solveit_evp/solveit_evp_mod1.f90
    libfepx/solveit_evp/solveit_evp_mod2.f90
    libfepx/solveit_evp/solveit_evp_mod3.f90
    libfepx/solveit_evp/solveit_evp_mod4.f90

    libfepx/solveit_evp/solve_evp/solve_evp_mod1.f90
    libfepx/solveit_evp/solve_evp/solve_evp_mod2.f90

    libfepx/solveit_evp/solve_evp/solve_evp_vpstress/solve_evp_vpstress_mod1.f90
    libfepx/solveit_evp/solve_evp/solve_evp_vpstress/solve_evp_vpstress_mod2.f90
    libfepx/solveit_evp/solve_evp/solve_evp_vpstress/solve_evp_vpstress_mod3.f90

    libfepx/solveit_evp/solve_evp/solve_evp_stress/solve_evp_stress_mod1.f90
    libfepx/solveit_evp/solve_evp/solve_evp_stress/solve_evp_stress_mod2.f90
    libfepx/solveit_evp/solve_evp/solve_evp_stress/solve_evp_stress_mod3.f90

    libfepx/solveit_evp/solve_evp/solve_evp_rstar/solve_evp_rstar_mod1.f90
    libfepx/solveit_evp/solve_evp/solve_evp_rstar/solve_evp_rstar_mod2.f90
    libfepx/solveit_evp/solve_evp/solve_evp_rstar/update_hardening/hardening_mod1.f90
    libfepx/solveit_evp/solve_evp/solve_evp_rstar/update_hardening/hardening_mod2.f90
    libfepx/solveit_evp/solve_evp/solve_evp_rstar/update_hardening/hardening_mod3.f90
    libfepx/solveit_evp/solve_evp/solve_evp_rstar/update_hardening/hardening_mod4.f90

    libfepx/solveit_evp/finalize_res/finalize_res_mod1.f90

    libfepx/restart/restart.f90

    libfepx/aniso_vp_lag.f90
    libfepx/kinematics_mod.f90
    libfepx/kinematics_mod_bis.f90
    
    libfepx/conjugate_gradient/conjugate_gradient_mod.f90
    libfepx/conjugate_gradient/conjugate_gradient_mod2.f90

    general/matrix_operations_mod.f90
    general/array_operations_mod.f90
    general/shape_2d_mod.f90
    general/shape_3d_mod.f90
    general/quadrature_mod.f90
    general/units_mod.f90

    libfepx/stiffness_mod.f90
    libfepx/surf_mod.f90
    libfepx/write_res/write_res_mod1.f90
    libfepx/write_res/write_res_mod2.f90
)

include_directories(${CMAKE_Fortran_MODULE_DIRECTORY} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

# Prepare MPI Fortan compiler for compilation of executable
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
message(STATUS "Switching to MPI Fortran compiler: ${CMAKE_Fortran_COMPILER}, ${CMAKE_Fortran_COMPILER_ID}")

message(STATUS "Preparing executable: fepx")
add_executable(fepx ${fepx})

# Install FEPX into bin on the local machine
install(TARGETS fepx DESTINATION bin
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
    GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Notify the user how to continue after `cmake ..`
message(STATUS "FEPX is set to be installed in: ${CMAKE_INSTALL_PREFIX}/bin")

# Testing ######################################################################

include(CTest)

file(GLOB_RECURSE tests ../tests/*/*.cmake)

if(BUILD_TESTING)
  foreach(test ${tests})
    get_filename_component(test_name ${test} NAME)
    string(REPLACE ".cmake" "" test_name ${test_name})
    get_filename_component(test_dir ${test} DIRECTORY)
    get_filename_component(test_path ${test} PATH)
    string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/../tests/" "" test_label ${test_dir})
    add_test(NAME ${test_label}
             COMMAND ${CMAKE_COMMAND}
             -Dtest_prog=$<TARGET_FILE:fepx>
             -Dtest_mode=${BUILD_TESTING_MODE}
             -Dtest_diff=${BUILD_TESTING_DIFF}
             -P ${test}
             WORKING_DIRECTORY ${test_dir})
  endforeach()
endif()
