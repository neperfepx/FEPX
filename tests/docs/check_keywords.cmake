# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Check that a generated page of the manual is the one its generator writes from the source:
# the slip-system lists (doc/slip_systems_table.rst, fepx_slip_systems, step 6.1; the
# keyword pages: check_keyword_pages.cmake). On failure, regenerate it with
#   <build>/src/<generator> <committed file>
#
# Variables: GENERATOR (the generator executable), COMMITTED (the committed page), GENERATED
# (scratch file in the build tree).

# Generate the reference from the current source table.
execute_process(COMMAND ${GENERATOR} ${GENERATED} RESULT_VARIABLE status)
if (NOT status EQUAL 0)
  message(FATAL_ERROR "${GENERATOR} failed (status ${status})")
endif()

# Compare with the committed copy.
execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files ${GENERATED} ${COMMITTED}
                RESULT_VARIABLE differ)
if (NOT differ EQUAL 0)
  message(FATAL_ERROR "${COMMITTED} is out of date with the source; regenerate it with: "
                      "${GENERATOR} ${COMMITTED}")
endif()
