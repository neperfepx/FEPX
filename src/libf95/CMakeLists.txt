# This file is part of the FEPX software package.
# Copyright (C) 1996-2020, DPLab, ACME Lab.
# See the COPYING file in the top-level directory.

message(STATUS "Preparing library: libf95")

file(GLOB libf95_src
    command_mod.f90
    files_mod.f90
    intrinsic_types_mod.f90
    libf95.f90
    strings_mod.f90
    timer_mod.f90
    )

add_library(libf95 STATIC ${libf95_src})
set_target_properties(libf95 PROPERTIES OUTPUT_NAME "f95dplab")
set_target_properties(libf95 PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(libf95 PROPERTIES Fortran_FORMAT "FREE")

install(TARGETS libf95 DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
