# This file is part of the FEPX software package.
# Copyright (C) 1996-2020, DPLab, ACME Lab.
# See the COPYING file in the top-level directory.

message(STATUS "Preparing library: libparallel")

file(GLOB libparallel_src
    gather_scatter_mod.f90
    parallel_mod.f90
    )

add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH} ${CMAKE_Fortran_MODULE_DIRECTORY})
link_directories(${MPI_Fortran_LIBRARIES})

add_library(libparallel STATIC ${libparallel_src})

add_dependencies(libparallel libf95)
set_target_properties(libparallel PROPERTIES OUTPUT_NAME "parallel")
set_target_properties(libparallel PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
set_target_properties(libparallel PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(libparallel PROPERTIES Fortran_FORMAT "FREE" )

install(TARGETS libparallel DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
