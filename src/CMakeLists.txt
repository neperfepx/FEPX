# This file is part of the FEPX software package.
# Copyright (C) 1996-2021, DPLab, ACME Lab.
# See the COPYING file in the top-level directory.

# Top-level CMakeLists.txt for building FEPX
cmake_minimum_required(VERSION 3.0)
project(fepx Fortran)
set(FEPX_VERSION \"1.2.2-1\")

# Prepare CMake for Fortran and locate MPI package
enable_language(Fortran)
find_package(MPI REQUIRED)

# Set the build type for appropriate compile flags
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build: Debug Release" FORCE)
endif()

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
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fcheck=all -fbacktrace -fbounds-check -Wall")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ftree-vectorize -funroll-loops")
# IFORT is currently not supported, but will be in the future
#ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
#   set(CMAKE_Fortran_MODDIR_FLAG   "")
#   set(CMAKE_Fortran_FLAGS         "")
#   set(CMAKE_Fortran_FLAGS_DEBUG   "")
#   set(CMAKE_Fortran_FLAGS_RELEASE "")
endif()

# Set directories for CMake to access files
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# Preamble info
STRING(TIMESTAMP CONFIG_DATE "%Y-%m-%d")
message(STATUS "Configuring Makefiles for the FEPX software package...")
message(STATUS "Built date: (${CONFIG_DATE})")
message(STATUS "OS: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}") # Need to grab input cache here?
message(STATUS "Version: ${VERSION}")

# Add in the subdirectories and includes
add_subdirectory(libf95)
add_subdirectory(libparallel)
add_subdirectory(libfepx)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

# Prepare MPI Fortan compiler for compilation of executable
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
message(STATUS "Switching to MPI Fortran compiler: ${CMAKE_Fortran_COMPILER}, ${CMAKE_Fortran_COMPILER_ID}")

message(STATUS "Preparing executable: fepx")
add_executable(fepx fepx.f90)
target_link_libraries(fepx PUBLIC libfepx libf95 libparallel)

# Install FEPX into bin on the local machine
install(TARGETS fepx DESTINATION bin
    PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ
    GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Notify the user how to continue after `cmake ..`
message(STATUS "FEPX is set to be installed in: ${CMAKE_INSTALL_PREFIX}/bin")
