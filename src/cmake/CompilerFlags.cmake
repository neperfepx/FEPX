# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Compiler flag sets for FEPX (conventions.md §6).
#
# Supported compilers:
#   GNU       gfortran (tested in CI)
#   IntelLLVM ifx      (flags provided; testing deferred, ROADMAP D6)
#
# All configurations: Fortran 2018 standard checking, implicit none, full warnings.
# Debug adds run-time checks, backtraces, floating-point traps (invalid, divide by zero,
# overflow) and signalling-NaN initialization of reals, so that uninitialized values fail loudly,
# at -O0. Checked adds the same run-time checks at -O2 (CMAKE_Fortran_FLAGS_CHECKED, top-level
# CMakeLists.txt; ROADMAP D88).
# Release keeps CMake's default optimization (-O3) and adds no value-changing flags
# (no -ffast-math, no -march=native), so that results are reproducible across machines.
# On arm64 (Apple silicon, aarch64 Linux), GNU builds also get -ffp-contract=off: there the base
# instruction set has fused multiply-add, which gfortran uses by default (-ffp-contract=fast), so
# results would differ from x86-64 in the last bits. x86-64 builds are unchanged (no FMA without
# -march; ROADMAP D112).

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

  # Standard conformance and warnings (all configurations).
  add_compile_options(
    -std=f2018
    -fimplicit-none
    -pedantic
    -Wall
    -Wextra
    -Wimplicit-interface
    -Wimplicit-procedure
  )
  # No optimization for Debug builds; run-time checks for Debug and Checked builds.
  add_compile_options(
    "$<$<CONFIG:Debug>:-O0;-g>"
    "$<$<CONFIG:Debug,Checked>:-fcheck=all;-fbacktrace;-ffpe-trap=invalid,zero,overflow;-finit-real=snan>"
  )
  # No fused multiply-add contraction on arm64 (results as on x86-64; see the header).
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64|ARM64)$")
    add_compile_options(-ffp-contract=off)
  endif()
  if(FEPX_WERROR)
    add_compile_options(-Werror)
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")

  # Standard conformance and warnings (all configurations).
  # (SHELL: keeps each option and its argument together; CMake would otherwise de-duplicate
  #  repeated words such as "all".)
  add_compile_options(
    "SHELL:-stand f18"
    "SHELL:-warn all"
  )
  # No optimization for Debug builds; run-time checks for Debug and Checked builds (argument
  # temporaries are allowed: they are not errors).
  add_compile_options(
    "$<$<CONFIG:Debug>:-O0;-g>"
    "$<$<CONFIG:Debug,Checked>:-traceback;-fpe0;-init=snan,arrays>"
    "$<$<CONFIG:Debug,Checked>:SHELL:-check all,noarg_temp_created>"
  )
  if(FEPX_WERROR)
    add_compile_options("SHELL:-warn errors")
  endif()

else()

  message(WARNING
    "Fortran compiler ${CMAKE_Fortran_COMPILER_ID} is not one of the supported compilers "
    "(GNU gfortran, IntelLLVM ifx). No FEPX-specific flags are set.")

endif()
