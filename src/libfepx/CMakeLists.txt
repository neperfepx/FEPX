# This file is part of the FEPX software package.
# Copyright (C) 1996-2020, DPLab, ACME Lab.
# See the COPYING file in the top-level directory.

message(STATUS "Preparing library: libfepx")

include_directories(${CMAKE_Fortran_MODULE_DIRECTORY} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

file(GLOB libfepx_src
    aniso_evps_mod.f90
    aniso_vp_lag.f90
    boundary_conditions_mod.f90
    conjugate_gradient_mod.f90
    convergence_mod.f90
    crystal_type_mod.f90
    dimensions_mod.f90
    driver_iso_vp_mod.f90
    driver_triaxclr_mod.f90
    driver_triaxcsr_mod.f90
    driver_uniaxial_control_mod.f90
    driver_utilities_mod.f90
    fiber_average_mod.f90
    hardening_mod.f90
    iterate_stress_evps_mod.f90
    iterate_stress_vp_mod.f90
    kinematics_mod.f90
    material_matrix_evps_mod.f90
    material_matrix_vp_mod.f90
    matrix_operations_mod.f90
    microstructure_mod.f90
    orientation_conversion_mod.f90
    polycrystal_response_evps_mod.f90
    quadrature_mod.f90
    read_input_mod.f90
    rstarn_solve_mod.f90
    shape_2d_mod.f90
    shape_3d_mod.f90
    stiffness_evps_mod.f90
    stiffness_vp_mod.f90
    stress_solve_evps_mod.f90
    stress_solve_vp_mod.f90
    surface_mod.f90
    units_mod.f90
    write_output_mod.f90
    )

add_library(libfepx STATIC ${libfepx_src})

add_dependencies(libfepx libf95 libparallel)
set_target_properties(libfepx PROPERTIES OUTPUT_NAME "fepx")
set_target_properties(libfepx PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(libfepx PROPERTIES Fortran_FORMAT "FREE")

install(TARGETS libfepx DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
