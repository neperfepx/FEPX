# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Run a command and check its exit status and output (CTest script, run with cmake -P).
#
# Variables (-D):
#   COMMAND          the command to run (a CMake list: program and arguments)
#   EXPECTED_STATUS  required exit status (integer), or "nonzero" for any failing status (used
#                    for runs through an MPI launcher ended by MPI_Abort: the launcher's exit
#                    status is implementation-defined)
#   EXPECTED_STDOUT  regular expression that standard output must match (optional; not
#                    checked when undefined or empty)
#   EXPECTED_STDERR  regular expression that standard error must match (optional, as above)
#   EXPECTED_STDOUT_ONCE  regular expression that standard output must match exactly once
#                    (optional; e.g. an error message written by rank 0 only, step 9.7)
#   WORKING_DIR      working directory of the command (optional)
#   REMOVE_BEFORE    a file or directory removed before the run (optional; e.g. the output
#                    directory of a previous run, which a new run refuses to overwrite)
#   EXPECTED_FILE    a file that must exist after the run, relative to WORKING_DIR (optional;
#                    e.g. the step results written before a run stops with an error)
#
# Used for the error_mod tests, which must check both the message and the exact exit status
# (PASS_REGULAR_EXPRESSION alone ignores the exit status, WILL_FAIL accepts any failure).
#
# The log's messages are wrapped at 72 columns (step 9.7, log_mod): before the regular
# expressions are matched, every continuation line is joined to its message with one blank
# (fepx_unwrap_log), and every run of blanks is made one blank, in the output and in the
# expressions alike (a break where a message had two blanks leaves one), so that a message
# matches whatever its line breaks.

include(${CMAKE_CURRENT_LIST_DIR}/unwrap_log.cmake)

if(NOT DEFINED WORKING_DIR)
  set(WORKING_DIR ${CMAKE_CURRENT_BINARY_DIR})
endif()
if(DEFINED REMOVE_BEFORE AND NOT REMOVE_BEFORE STREQUAL "")
  file(REMOVE_RECURSE ${REMOVE_BEFORE})
endif()

execute_process(COMMAND ${COMMAND}
                WORKING_DIRECTORY ${WORKING_DIR}
                RESULT_VARIABLE status
                OUTPUT_VARIABLE stdout
                ERROR_VARIABLE stderr)

message("standard output:\n${stdout}\nstandard error:\n${stderr}")

# The messages joined across their line breaks (wrapping at 72 columns, step 9.7).
fepx_unwrap_log("${stdout}" stdout)
fepx_unwrap_log("${stderr}" stderr)
string(REGEX REPLACE "  +" " " stdout "${stdout}")
string(REGEX REPLACE "  +" " " stderr "${stderr}")
string(REGEX REPLACE "  +" " " EXPECTED_STDOUT "${EXPECTED_STDOUT}")
string(REGEX REPLACE "  +" " " EXPECTED_STDERR "${EXPECTED_STDERR}")

if(EXPECTED_STATUS STREQUAL "nonzero")
  if(status STREQUAL "0" OR NOT status MATCHES "^[0-9]+$")
    message(FATAL_ERROR "exit status: actual '${status}', expected a non-zero status")
  endif()
elseif(NOT status STREQUAL "${EXPECTED_STATUS}")
  message(FATAL_ERROR "exit status: actual '${status}', expected '${EXPECTED_STATUS}'")
endif()

if(NOT "${EXPECTED_STDOUT}" STREQUAL "" AND NOT stdout MATCHES "${EXPECTED_STDOUT}")
  message(FATAL_ERROR "standard output does not match the regular expression "
                      "'${EXPECTED_STDOUT}'")
endif()

if(NOT "${EXPECTED_STDOUT_ONCE}" STREQUAL "")
  string(REGEX REPLACE "  +" " " EXPECTED_STDOUT_ONCE "${EXPECTED_STDOUT_ONCE}")
  string(REGEX MATCHALL "${EXPECTED_STDOUT_ONCE}" matches "${stdout}")
  list(LENGTH matches num_matches)
  if(NOT num_matches EQUAL 1)
    message(FATAL_ERROR "standard output matches the regular expression "
                        "'${EXPECTED_STDOUT_ONCE}' ${num_matches} times, not once")
  endif()
endif()

if(NOT "${EXPECTED_STDERR}" STREQUAL "" AND NOT stderr MATCHES "${EXPECTED_STDERR}")
  message(FATAL_ERROR "standard error does not match the regular expression "
                      "'${EXPECTED_STDERR}'")
endif()

if(DEFINED EXPECTED_FILE AND NOT EXPECTED_FILE STREQUAL "" AND
   NOT EXISTS ${WORKING_DIR}/${EXPECTED_FILE})
  message(FATAL_ERROR "the run did not write ${EXPECTED_FILE}")
endif()
