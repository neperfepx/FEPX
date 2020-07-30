# This file is part of the FEPX software package.
# Copyright (C) 1996-2020, DPLab, ACME Lab.
# See the COPYING file in the top-level directory.

message(STATUS "Preparing library: libfepx")

include_directories(${CMAKE_Fortran_MODULE_DIRECTORY} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

file(GLOB libfepx_src
    aniso_evps_lag.f90
    aniso_vp_lag.f90
    boundary_conditions_mod.f90
    ConvergenceModule.f90
    CrystalTypeModule.f90
    DimsModule.f90
    Dpeff.f90
    DriverVpModule.f90
    driver_uniaxial_control_mod.f90
    driver_triaxclr_mod.f90
    driver_triaxcsr_mod.f90
    driver_utilities_mod.f90
    elemental_variables_utils_mod.f90
    F90UtilitiesModule.f90
    Hardening.f90
    hart.f90
    itstress_evps_lag.f90
    itstress_vp_lag.f90
    KinematicsModule.f90
    material_mtx_evps_lag.f90
    material_mtx_vp_lag.f90
    read_input_mod.f90
    microstructure.f90
    orientation_conversion_mod.f90
    parallel-matrix.f90
    pcg_solver.f90
    poly_response_evps_lag.f90
    post_update_n.f90
    powder_diffraction_mod.f90
    quadrature.f90
    RestartModule.f90
    rstarn_solve_lag.f90
    shape_2d.f90
    shape_3d.f90
    simulation_configuration_mod.f90
    StiffnessEvpsModule.f90
    stiffness_vp_lag.f90
    stress_deform_lag.f90
    stress_solve_evps_lag.f90
    stress_solve_vp_lag.f90
    stress_strain.f90
    surface.f90
    surf_info.f90
    Tensor3DModule.f90
    units.f90
    utils_crystal3D_lag.f90
    write_output_mod.f90
    )

add_library(libfepx STATIC ${libfepx_src})

add_dependencies(libfepx libf95 libparallel)
set_target_properties(libfepx PROPERTIES OUTPUT_NAME "fepx")
set_target_properties(libfepx PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(libfepx PROPERTIES Fortran_FORMAT "FREE")

install(TARGETS libfepx DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
