# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Installation of FEPX for the user, without administrator rights (CMake script, run with cmake -P
# by the target userinstall of src/cmake/Packaging.cmake; as Neper's make userinstall).
#
# Variables (-D):
#   BINARY_DIR  the build directory (configured and built)
#
# Installs the build, as cmake --install does, under $HOME/.local, the user's prefix of the XDG
# base directories (bin/fepx, share/man/man1/fepx.1.gz, share/doc/fepx/), whatever
# CMAKE_INSTALL_PREFIX is. HOME is read when the target runs (not when CMake configures), so that
# the installation goes to the home directory of the user who runs it. A warning says when
# $HOME/.local/bin is not in the PATH.
#
# Local variables:
#   home      the home directory ($HOME)
#   prefix    the installation prefix, $HOME/.local
#   result    the exit status of cmake --install
#   path_dirs the directories of the PATH, as a list

# The home directory, required.
set(home "$ENV{HOME}")
if(home STREQUAL "")
  message(FATAL_ERROR "userinstall: HOME is not set.")
endif()
set(prefix "${home}/.local")

# The installation, as cmake --install with the prefix $HOME/.local.
message(STATUS "userinstall: installing into ${prefix}")
execute_process(COMMAND "${CMAKE_COMMAND}" --install "${BINARY_DIR}" --prefix "${prefix}"
                RESULT_VARIABLE result)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "userinstall: cmake --install failed.")
endif()

# A warning when the executable's directory is not in the PATH.
if(CMAKE_HOST_WIN32)
  set(path_dirs "$ENV{PATH}")
else()
  string(REPLACE ":" ";" path_dirs "$ENV{PATH}")
endif()
if(NOT "${prefix}/bin" IN_LIST path_dirs AND NOT "${prefix}/bin/" IN_LIST path_dirs)
  message(WARNING "userinstall: ${prefix}/bin is not in the PATH; add it (e.g. in ~/.bashrc: "
                  "export PATH=\"$HOME/.local/bin:$PATH\") to run fepx from any directory.")
endif()
message(STATUS "userinstall: fepx installed as ${prefix}/bin/fepx")
