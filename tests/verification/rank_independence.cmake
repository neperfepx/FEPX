# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Run one v2-fixed case on several rank counts and check that the results agree
# (CLAUDE.md §6; CHANGES.md C36), CMake script mode (`cmake -P`).
#
# Each run gets a fresh directory WORK_DIR/np<n> with the case inputs (CASE_DIR/simulation.cfg,
# simulation.opt when present) and MESH as simulation.msh; then CHECK_SCRIPT compares the element stresses
# of every run with those of the first (relative tolerance RTOL).
#
# Variables (passed with -D): EXECUTABLE, CASE_DIR, MESH, WORK_DIR, MPIEXEC, MPIEXEC_NUMPROC,
# RANKS (list, "|" separated), PYTHON, CHECK_SCRIPT, RTOL, MPIEXEC_EXTRA (space separated,
# optional), RUN_TIMEOUT (time limit of each run [s], optional, default 240).

foreach(var EXECUTABLE CASE_DIR MESH WORK_DIR MPIEXEC MPIEXEC_NUMPROC RANKS PYTHON CHECK_SCRIPT
            RTOL)
  if(NOT DEFINED ${var})
    message(FATAL_ERROR "v2_rank_independence.cmake: ${var} is not set")
  endif()
endforeach()
separate_arguments(extra_args UNIX_COMMAND "${MPIEXEC_EXTRA}")
if(NOT DEFINED RUN_TIMEOUT)
  set(RUN_TIMEOUT 240)
endif()
string(REPLACE "|" ";" rank_list "${RANKS}")

# One run per rank count.
set(run_dirs "")
foreach(n IN LISTS rank_list)
  set(dir ${WORK_DIR}/np${n})
  file(REMOVE_RECURSE ${dir})
  file(MAKE_DIRECTORY ${dir})
  configure_file(${CASE_DIR}/simulation.cfg ${dir}/simulation.cfg COPYONLY)
  if(EXISTS ${CASE_DIR}/simulation.opt)
    configure_file(${CASE_DIR}/simulation.opt ${dir}/simulation.opt COPYONLY)
  endif()
  # The mesh's copy with n partitions when it has one (FEPXTesting.cmake, fepx_mesh_for_ranks).
  set(mesh_for_ranks ${MESH})
  if(NOT n EQUAL 1)
    string(REGEX REPLACE "(_part[0-9]+)?[.]msh$" "" mesh_stem "${MESH}")
    foreach(candidate "${mesh_stem}_part${n}.msh" "${mesh_stem}_blocks${n}.msh")
      if(EXISTS ${candidate})
        set(mesh_for_ranks ${candidate})
        break()
      endif()
    endforeach()
  endif()
  configure_file(${mesh_for_ranks} ${dir}/simulation.msh COPYONLY)
  execute_process(COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC} ${n} ${extra_args} ${EXECUTABLE}
                  WORKING_DIRECTORY ${dir}
                  TIMEOUT ${RUN_TIMEOUT}
                  RESULT_VARIABLE status
                  OUTPUT_VARIABLE output
                  ERROR_VARIABLE output)
  if(NOT status STREQUAL "0")
    message(FATAL_ERROR "run on ${n} rank(s) failed (status '${status}'):\n${output}")
  endif()
  list(APPEND run_dirs ${dir})
endforeach()

# Compare every run with the first.
execute_process(COMMAND ${PYTHON} ${CHECK_SCRIPT} --rtol=${RTOL} ${run_dirs}
                RESULT_VARIABLE status)
if(NOT status STREQUAL "0")
  message(FATAL_ERROR "results depend on the number of ranks")
endif()
