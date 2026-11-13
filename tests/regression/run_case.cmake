# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Run one regression case and compare it with its reference data (CMake script mode, `cmake -P`).
#
# Source : FEPX regression harness (ROADMAP step 0.3); policy in tests/regression/README.md.
#
# Steps:
#   1. Recreate WORK_DIR and copy the case inputs into it (simulation.cfg, .msh, .opt, .tess,
#      whichever exist; MESH_FILE, ORI_FILE and PHASE_FILE, if given, are copied as
#      simulation.msh, .ori and .phase). The case source directory is never written to.
#   2. If PREREQ_SIM_DIR is set (second part of a restart pair), copy that simulation.sim
#      directory (its results and restart point) into WORK_DIR.
#   3. Run the executable with MPI on NUM_RANKS ranks; it must exit with status 0 and write
#      WORK_DIR/simulation.sim/results.
#   3b. If CHECK_SCRIPT is set (verification cases), run `PYTHON CHECK_SCRIPT WORK_DIR/simulation.sim
#      CHECK_ARGS`; it must exit with status 0.
#   4. Unless COMPARE is OFF:
#      - normal mode: compare WORK_DIR/simulation.sim with REFERENCE_DIR using
#        tests/tools/compare_sim.py (COMPARE_ARGS are passed through);
#      - update mode (environment variable FEPX_REGRESSION_UPDATE): regenerate the reference
#        instead of failing, under the rules of CLAUDE.md §2.5:
#          =1    replace the reference only if the new results fail the comparison (or there is
#                no reference), so that run-to-run noise below the tolerance (e.g. the last printed
#                digit on 4 ranks) does not churn the reference data;
#          =all  replace every reference unconditionally.
#
# Variables (passed with -D):
#   EXECUTABLE       program to run (required)
#   CASE_SOURCE_DIR  directory holding the case inputs (required)
#   WORK_DIR         run directory in the build tree, recreated on every run (required)
#   MPIEXEC          MPI launcher (required)
#   MPIEXEC_NUMPROC  launcher option giving the number of ranks, e.g. -n (required)
#   NUM_RANKS        number of MPI ranks (required)
#   MESH_FILE        mesh to use as simulation.msh (optional; default CASE_SOURCE_DIR's own)
#   ORI_FILE, PHASE_FILE  files to use as simulation.ori, simulation.phase (optional)
#   MPIEXEC_EXTRA    extra launcher options, space separated (optional)
#   PREREQ_SIM_DIR   simulation.sim directory of the run a restart continues (optional)
#   COMPARE          ON (default) or OFF (run only)
#   REFERENCE_DIR    reference .sim directory (required if COMPARE is ON)
#   REFERENCE_READ_ONLY  if true, update mode never writes REFERENCE_DIR (the reference belongs to
#                    another case; this case only compares)
#   PYTHON           Python 3 interpreter (required if COMPARE is ON)
#   COMPARE_SCRIPT   path of tests/tools/compare_sim.py (required if COMPARE is ON)
#   COMPARE_ARGS     extra compare_sim.py options, separated by '|' (optional)
#   CHECK_SCRIPT     Python check script run after the simulation (optional; needs PYTHON)
#   CHECK_ARGS       its options after the .sim directory, separated by '|' (optional)

# --- Arguments ----------------------------------------------------------------------------------
foreach(var EXECUTABLE CASE_SOURCE_DIR WORK_DIR MPIEXEC MPIEXEC_NUMPROC NUM_RANKS)
  if(NOT DEFINED ${var})
    message(FATAL_ERROR "run_case.cmake: ${var} is not set")
  endif()
endforeach()
if(NOT DEFINED COMPARE)
  set(COMPARE ON)
endif()
if(COMPARE)
  foreach(var REFERENCE_DIR PYTHON COMPARE_SCRIPT)
    if(NOT DEFINED ${var})
      message(FATAL_ERROR "run_case.cmake: ${var} is not set (needed with COMPARE ON)")
    endif()
  endforeach()
endif()

# --- 1. Fresh copy of the inputs ----------------------------------------------------------------
file(REMOVE_RECURSE ${WORK_DIR})
file(MAKE_DIRECTORY ${WORK_DIR})
set(num_inputs 0)
foreach(input simulation.cfg simulation.msh simulation.opt simulation.tess)
  if(EXISTS ${CASE_SOURCE_DIR}/${input})
    file(COPY ${CASE_SOURCE_DIR}/${input} DESTINATION ${WORK_DIR})
    math(EXPR num_inputs "${num_inputs} + 1")
  endif()
endforeach()
foreach(input MESH_FILE:simulation.msh ORI_FILE:simulation.ori PHASE_FILE:simulation.phase)
  string(REPLACE ":" ";" input "${input}")
  list(GET input 0 input_var)
  list(GET input 1 input_name)
  if(DEFINED ${input_var} AND NOT ${input_var} STREQUAL "")
    file(COPY_FILE ${${input_var}} ${WORK_DIR}/${input_name})
  endif()
endforeach()
if(NOT EXISTS ${WORK_DIR}/simulation.cfg OR NOT EXISTS ${WORK_DIR}/simulation.msh)
  message(FATAL_ERROR "case ${CASE_SOURCE_DIR}: simulation.cfg or simulation.msh missing")
endif()

# --- 2. The simulation of the prerequisite case -------------------------------------------------
if(DEFINED PREREQ_SIM_DIR AND NOT PREREQ_SIM_DIR STREQUAL "")
  if(NOT IS_DIRECTORY ${PREREQ_SIM_DIR})
    message(FATAL_ERROR "no ${PREREQ_SIM_DIR} (run the prerequisite first)")
  endif()
  file(COPY ${PREREQ_SIM_DIR} DESTINATION ${WORK_DIR})
endif()

# --- 3. Run ---------------------------------------------------------------------------------------
separate_arguments(extra_args UNIX_COMMAND "${MPIEXEC_EXTRA}")
execute_process(
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC} ${NUM_RANKS} ${extra_args} ${EXECUTABLE}
  WORKING_DIRECTORY ${WORK_DIR}
  RESULT_VARIABLE run_status
  OUTPUT_FILE ${WORK_DIR}/run.log
  ERROR_FILE ${WORK_DIR}/run.log)
if(NOT run_status EQUAL 0)
  file(READ ${WORK_DIR}/run.log run_log)
  message(FATAL_ERROR "run failed (status ${run_status}):\n${run_log}")
endif()
# (The text results, or with output_format hdf5 alone the HDF5 file, step 9.1.)
if(NOT EXISTS ${WORK_DIR}/simulation.sim/results AND NOT EXISTS ${WORK_DIR}/simulation.sim/simulation.h5)
  message(FATAL_ERROR "the run did not write ${WORK_DIR}/simulation.sim/results")
endif()
message(STATUS "run completed: ${WORK_DIR}/simulation.sim")

# --- 3b. Verification check ------------------------------------------------------------------------
if(DEFINED CHECK_SCRIPT AND NOT CHECK_SCRIPT STREQUAL "")
  string(REPLACE "|" ";" check_args "${CHECK_ARGS}")
  execute_process(
    COMMAND ${PYTHON} ${CHECK_SCRIPT} ${WORK_DIR}/simulation.sim ${check_args}
    RESULT_VARIABLE check_status
    OUTPUT_VARIABLE check_output
    ERROR_VARIABLE check_output)
  message("${check_output}")
  if(NOT check_status EQUAL 0)
    message(FATAL_ERROR "verification check failed (${CHECK_SCRIPT})")
  endif()
endif()

if(NOT COMPARE)
  return()
endif()

# --- 4a. Compare with the reference (if there is one) ------------------------------------------
set(compare_status 1)
if(IS_DIRECTORY ${REFERENCE_DIR}/results)
  string(REPLACE "|" ";" compare_args "${COMPARE_ARGS}")
  execute_process(
    COMMAND ${PYTHON} ${COMPARE_SCRIPT} ${WORK_DIR}/simulation.sim ${REFERENCE_DIR} --quiet
            ${compare_args}
    RESULT_VARIABLE compare_status
    OUTPUT_VARIABLE compare_output
    ERROR_VARIABLE compare_output)
  message("${compare_output}")
endif()

# --- 4b. Update mode: replace the reference results ---------------------------------------------
if(REFERENCE_READ_ONLY AND NOT "$ENV{FEPX_REGRESSION_UPDATE}" STREQUAL "")
  if(NOT compare_status EQUAL 0)
    message(FATAL_ERROR "results differ from the (read-only) reference ${REFERENCE_DIR}")
  endif()
  return()
endif()
if("$ENV{FEPX_REGRESSION_UPDATE}" STREQUAL "all" OR
   ("$ENV{FEPX_REGRESSION_UPDATE}" STREQUAL "1" AND NOT compare_status EQUAL 0))
  file(REMOVE_RECURSE ${REFERENCE_DIR}/results)
  file(MAKE_DIRECTORY ${REFERENCE_DIR})
  file(COPY ${WORK_DIR}/simulation.sim/results DESTINATION ${REFERENCE_DIR})
  message(STATUS "reference updated: ${REFERENCE_DIR}")
  return()
elseif("$ENV{FEPX_REGRESSION_UPDATE}" STREQUAL "1")
  message(STATUS "reference unchanged (agrees within tolerance): ${REFERENCE_DIR}")
  return()
endif()

# --- 4c. Normal mode: the comparison decides -------------------------------------------------------
if(NOT IS_DIRECTORY ${REFERENCE_DIR}/results)
  message(FATAL_ERROR "no reference data in ${REFERENCE_DIR} "
                      "(generate it as described in tests/regression/README.md)")
endif()
if(NOT compare_status EQUAL 0)
  message(FATAL_ERROR "results differ from the reference ${REFERENCE_DIR}")
endif()
