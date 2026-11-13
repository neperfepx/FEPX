# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Restart errors of the fepx executable (step 8.4, dev/theory/addendum_restart.md §3, §7),
# CMake script mode (`cmake -P`): part 1 of a restart case is run once with print restart; each
# case then continues its output directory with a modified part 2 and must end with the given
# exit status and message:
#   ranks          part 2 on 1 rank (the restart files were written on 2)
#   material       tau_0 changed (a strict fingerprint line)
#   output_format  output_format ascii hdf5 (part 1: the default, ascii; a change of the output
#                  formats is an error, step 9.1; without HDF5 the configuration's error)
#   output_format_same  output_format ascii given (part 1: the default): the same output, runs
#                  to the end without a warning naming output_format
#   no_files       restart on without simulation.sim
#   no_restart_on  part 2 without restart on (the output directory exists)
#   warning        cg_tol changed: runs to the end, with a warning naming the line
#
# Variables (passed with -D): EXECUTABLE (fepx), CASE_DIR (part1.cfg, part2.cfg), MESH,
# WORK_DIR (removed and recreated), MPIEXEC, MPIEXEC_NUMPROC, MPIEXEC_EXTRA (optional), HDF5
# (ON or OFF: whether fepx has the HDF5 output).

foreach(var EXECUTABLE CASE_DIR MESH WORK_DIR MPIEXEC MPIEXEC_NUMPROC)
  if(NOT DEFINED ${var})
    message(FATAL_ERROR "restart_errors.cmake: ${var} is not set")
  endif()
endforeach()
separate_arguments(extra_args UNIX_COMMAND "${MPIEXEC_EXTRA}")

# The joining of the log's wrapped messages (step 9.7).
include(${CMAKE_CURRENT_LIST_DIR}/../core/unwrap_log.cmake)

# Part 1, once.
file(REMOVE_RECURSE ${WORK_DIR})
file(MAKE_DIRECTORY ${WORK_DIR}/base)
configure_file(${MESH} ${WORK_DIR}/base/simulation.msh COPYONLY)
configure_file(${CASE_DIR}/part1.cfg ${WORK_DIR}/base/simulation.cfg COPYONLY)
execute_process(COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC} 2 ${extra_args} ${EXECUTABLE}
                WORKING_DIRECTORY ${WORK_DIR}/base TIMEOUT 300 RESULT_VARIABLE status
                OUTPUT_VARIABLE output ERROR_VARIABLE output)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "part 1 failed (status ${status}):\n${output}")
endif()
file(READ ${CASE_DIR}/part2.cfg part2)

set(num_failed 0)
foreach(case ranks material output_format output_format_same no_files no_restart_on warning)
  # The case's directory: the mesh, part 1's output (except no_files), the modified part 2.
  set(dir ${WORK_DIR}/${case})
  file(MAKE_DIRECTORY ${dir})
  configure_file(${MESH} ${dir}/simulation.msh COPYONLY)
  if(NOT case STREQUAL "no_files")
    file(COPY ${WORK_DIR}/base/simulation.sim DESTINATION ${dir})
  endif()
  set(cfg "${part2}")
  set(ranks 2)
  set(expect_ok FALSE)
  set(absent "")
  if(case STREQUAL "ranks")
    set(ranks 1)
    set(pattern "written on 2 rank[(]s[)]. a restart needs the same number of ranks")
  elseif(case STREQUAL "material")
    string(REPLACE "tau_0 210.0" "tau_0 200.0" cfg "${cfg}")
    set(pattern "\"tau_0 phase 1: 210.0\" of the restart point's configuration is not in simulation.cfg")
  elseif(case STREQUAL "output_format")
    set(cfg "${cfg}output_format hdf5 ascii\n")
    if(HDF5)
      set(pattern "\"output_format: ascii\" of the restart point's configuration is not in simulation.cfg [(]the material, the boundary conditions, the input files and the output format must not change[)]")
    else()
      set(pattern "output_format hdf5: this FEPX build has no HDF5 output")
    endif()
  elseif(case STREQUAL "output_format_same")
    set(cfg "${cfg}output_format ascii\n")
    set(pattern "Simulation completed successfully")
    set(absent "output_format")
    set(expect_ok TRUE)
  elseif(case STREQUAL "no_files")
    set(pattern "restart on: no restart file simulation.sim/restart/restart.rank0")
  elseif(case STREQUAL "no_restart_on")
    string(REPLACE "restart on\n" "" cfg "${cfg}")
    set(pattern "output directory simulation.sim already exists")
  else()
    set(cfg "${cfg}cg_tol 1e-7\n")
    set(pattern "\"cg_tol: 1e-7\" in simulation.cfg is not in the restart point's configuration")
    set(expect_ok TRUE)
  endif()
  file(WRITE ${dir}/simulation.cfg "${cfg}")

  # Run, then the status and the message.
  execute_process(COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC} ${ranks} ${extra_args} ${EXECUTABLE}
                  WORKING_DIRECTORY ${dir} TIMEOUT 300 RESULT_VARIABLE status
                  OUTPUT_VARIABLE output ERROR_VARIABLE output)
  # (The messages joined across their line breaks, blanks collapsed: unwrap_log.cmake.)
  fepx_normalize_log("${output}" output)
  string(REGEX MATCH "${pattern}" found "${output}")
  if(NOT absent STREQUAL "")
    string(FIND "${output}" "${absent}" position)
    if(NOT position EQUAL -1)
      set(found "")
    endif()
  endif()
  if(expect_ok)
    set(status_ok FALSE)
    if(status EQUAL 0)
      set(status_ok TRUE)
    endif()
  else()
    set(status_ok FALSE)
    if(status MATCHES "^[0-9]+$" AND NOT status EQUAL 0)
      set(status_ok TRUE)
    endif()
  endif()
  if(status_ok AND found)
    message(STATUS "ok    ${case}: exit status ${status}, message found")
  else()
    message(STATUS "FAIL  ${case}: exit status '${status}', message found: '${found}'\n${output}")
    math(EXPR num_failed "${num_failed} + 1")
  endif()
endforeach()

if(num_failed GREATER 0)
  message(FATAL_ERROR "${num_failed} restart error check(s) failed")
endif()
