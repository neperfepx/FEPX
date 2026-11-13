# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# The terminal-output golden-file test (step 9.7, ROADMAP D94; dev/plan/terminal_output.md §3.8;
# CMake script mode, cmake -P): run the case of CASE_DIR (simulation.cfg, with MESH) through the
# launcher command COMMAND in WORK_DIR, then compare its standard output with REFERENCE after
# masking the variable fields (CHECK_SCRIPT, check_terminal_output.py). The run must end with
# exit status 0; its standard error (the launcher's own lines) is not checked.
#
# Variables (-D): COMMAND (a CMake list: the launcher, its arguments and fepx), CASE_DIR, MESH,
# WORK_DIR (removed and recreated), PYTHON, CHECK_SCRIPT, REFERENCE.

foreach(var COMMAND CASE_DIR MESH WORK_DIR PYTHON CHECK_SCRIPT REFERENCE)
  if(NOT DEFINED ${var})
    message(FATAL_ERROR "terminal_output.cmake: ${var} is not set")
  endif()
endforeach()

# The case's directory.
file(REMOVE_RECURSE ${WORK_DIR})
file(MAKE_DIRECTORY ${WORK_DIR})
configure_file(${CASE_DIR}/simulation.cfg ${WORK_DIR}/simulation.cfg COPYONLY)
configure_file(${MESH} ${WORK_DIR}/simulation.msh COPYONLY)

# The run, its standard output kept in a file.
execute_process(COMMAND ${COMMAND}
                WORKING_DIRECTORY ${WORK_DIR}
                RESULT_VARIABLE status
                OUTPUT_FILE ${WORK_DIR}/stdout.txt
                ERROR_FILE ${WORK_DIR}/stderr.txt)
file(READ ${WORK_DIR}/stdout.txt stdout)
message("standard output:\n${stdout}")
if(NOT status EQUAL 0)
  file(READ ${WORK_DIR}/stderr.txt stderr)
  message(FATAL_ERROR "run failed (status ${status}):\n${stderr}")
endif()

# The comparison with the reference (masked).
execute_process(COMMAND ${PYTHON} ${CHECK_SCRIPT} ${WORK_DIR}/stdout.txt ${REFERENCE}
                RESULT_VARIABLE check_status)
if(NOT check_status EQUAL 0)
  message(FATAL_ERROR "the terminal output differs from ${REFERENCE}")
endif()
