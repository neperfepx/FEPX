# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Restart chain test (step 8.4, dev/theory/addendum_restart.md §8), CMake script mode
# (`cmake -P`): a run of the whole loading history (full.cfg) and the same history interrupted and
# continued from restart files (part1.cfg, part2.cfg, ... run one after the other in one
# directory; each part after the first has restart on and lists only its own steps) must leave
# byte-identical simulation.sim directories: every result file (steps, forces, convergence,
# macroscopic file) the same, apart from the copied inputs (inputs/) and the
# restart files (restart/), which hold the run's own configuration. The HDF5 file
# (output_format hdf5, step 9.1) is compared by h5diff (every dataset and attribute the same;
# its bytes differ: the datasets removed at a restart leave free space in the file).
#
# Variables (passed with -D): EXECUTABLE (fepx), CASE_DIR (full.cfg, part<i>.cfg, and optional
# auxiliary inputs simulation.ori/.phase/.opt), MESH (the mesh, copied as simulation.msh),
# WORK_DIR (removed and recreated), MPIEXEC, MPIEXEC_NUMPROC, MPIEXEC_EXTRA (optional),
# NUM_RANKS, RUN_TIMEOUT (seconds per run; optional, default 300), H5DIFF (h5diff; needed when
# the runs write simulation.h5).

foreach(var EXECUTABLE CASE_DIR MESH WORK_DIR MPIEXEC MPIEXEC_NUMPROC NUM_RANKS)
  if(NOT DEFINED ${var})
    message(FATAL_ERROR "restart_chain.cmake: ${var} is not set")
  endif()
endforeach()
if(NOT DEFINED RUN_TIMEOUT)
  set(RUN_TIMEOUT 300)
endif()
separate_arguments(extra_args UNIX_COMMAND "${MPIEXEC_EXTRA}")

# Run fepx in a directory; a failed run ends the test with its log.
function(run_fepx dir label)
  execute_process(COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC} ${NUM_RANKS} ${extra_args} ${EXECUTABLE}
                  WORKING_DIRECTORY ${dir}
                  TIMEOUT ${RUN_TIMEOUT}
                  RESULT_VARIABLE status
                  OUTPUT_FILE ${dir}/${label}.log
                  ERROR_FILE ${dir}/${label}.log)
  if(NOT status EQUAL 0)
    file(READ ${dir}/${label}.log log)
    message(FATAL_ERROR "${label}: run failed (status ${status}):\n${log}")
  endif()
endfunction()

# Fresh directories with the mesh and the auxiliary inputs.
file(REMOVE_RECURSE ${WORK_DIR})
foreach(sub full chain)
  file(MAKE_DIRECTORY ${WORK_DIR}/${sub})
  configure_file(${MESH} ${WORK_DIR}/${sub}/simulation.msh COPYONLY)
  foreach(aux simulation.ori simulation.phase simulation.opt)
    if(EXISTS ${CASE_DIR}/${aux})
      configure_file(${CASE_DIR}/${aux} ${WORK_DIR}/${sub}/${aux} COPYONLY)
    endif()
  endforeach()
endforeach()

# The whole history.
configure_file(${CASE_DIR}/full.cfg ${WORK_DIR}/full/simulation.cfg COPYONLY)
run_fepx(${WORK_DIR}/full full)

# The parts, in order, in one directory; each part after the first must continue a restart point.
file(GLOB parts RELATIVE ${CASE_DIR} ${CASE_DIR}/part*.cfg)
list(SORT parts)
list(LENGTH parts num_parts)
if(num_parts LESS 2)
  message(FATAL_ERROR "${CASE_DIR}: at least part1.cfg and part2.cfg are needed")
endif()
set(i_part 0)
foreach(part ${parts})
  math(EXPR i_part "${i_part} + 1")
  configure_file(${CASE_DIR}/${part} ${WORK_DIR}/chain/simulation.cfg COPYONLY)
  run_fepx(${WORK_DIR}/chain part${i_part})
  if(i_part GREATER 1)
    file(READ ${WORK_DIR}/chain/part${i_part}.log log)
    string(FIND "${log}" "Restart: continuing from step" pos)
    if(pos EQUAL -1)
      message(FATAL_ERROR "part${i_part}: the run did not continue from a restart point")
    endif()
  endif()
endforeach()

# The two output directories: the same files, byte for byte (inputs/ and restart/ excluded).
foreach(sub full chain)
  file(GLOB_RECURSE files_${sub} RELATIVE ${WORK_DIR}/${sub}/simulation.sim
       ${WORK_DIR}/${sub}/simulation.sim/*)
  list(FILTER files_${sub} EXCLUDE REGEX "^(inputs|restart)/")
  list(SORT files_${sub})
endforeach()
if(NOT files_full STREQUAL files_chain)
  message(FATAL_ERROR "the output files differ:\n  full:  ${files_full}\n  chain: ${files_chain}")
endif()
set(num_differ 0)
foreach(f ${files_full})
  if(f STREQUAL "simulation.h5")
    if(NOT DEFINED H5DIFF OR H5DIFF STREQUAL "")
      message(FATAL_ERROR "the runs wrote simulation.h5: H5DIFF is needed")
    endif()
    execute_process(COMMAND ${H5DIFF} ${WORK_DIR}/full/simulation.sim/${f}
                            ${WORK_DIR}/chain/simulation.sim/${f}
                    RESULT_VARIABLE differs OUTPUT_VARIABLE h5diff_out ERROR_VARIABLE h5diff_out)
    if(NOT differs EQUAL 0)
      message(STATUS "h5diff:\n${h5diff_out}")
    endif()
  else()
    execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files ${WORK_DIR}/full/simulation.sim/${f}
                            ${WORK_DIR}/chain/simulation.sim/${f}
                    RESULT_VARIABLE differs)
  endif()
  if(NOT differs EQUAL 0)
    message(STATUS "differs: ${f}")
    math(EXPR num_differ "${num_differ} + 1")
  endif()
endforeach()
list(LENGTH files_full num_files)
if(num_differ GREATER 0)
  message(FATAL_ERROR "${num_differ} of ${num_files} output files differ between the whole run "
                      "and the restarted chain")
endif()
message(STATUS "${num_files} output files identical (${num_parts} parts)")
