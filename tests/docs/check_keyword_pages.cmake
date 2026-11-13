# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Check that the manual's generated keyword pages are those fepx_keywords writes from the source
# (task M.1, step 9.9, ROADMAP D104): the keyword reference (keywords_table.rst), the
# compatibility table (compatibility_table.rst) and the entry of every implemented keyword
# (generated/keywords/<name>.inc), with no committed entry left over. On failure, regenerate
# them with
#   <build>/fepx_keywords doc
#
# Variables: GENERATOR (fepx_keywords executable), MANUAL_DIR (the committed manual,
# doc), WORK_DIR (scratch directory in the build tree).

# Generate the pages into the scratch directory.
file(REMOVE_RECURSE ${WORK_DIR})
file(MAKE_DIRECTORY ${WORK_DIR}/generated/keywords)
execute_process(COMMAND ${GENERATOR} ${WORK_DIR} RESULT_VARIABLE status)
if (NOT status EQUAL 0)
  message(FATAL_ERROR "${GENERATOR} failed (status ${status})")
endif()

# The same entry files, generated and committed.
file(GLOB generated_entries RELATIVE ${WORK_DIR} ${WORK_DIR}/generated/keywords/*.inc)
file(GLOB committed_entries RELATIVE ${MANUAL_DIR} ${MANUAL_DIR}/generated/keywords/*.inc)
list(SORT generated_entries)
list(SORT committed_entries)
if (NOT generated_entries STREQUAL committed_entries)
  message(FATAL_ERROR "the committed keyword entries (${MANUAL_DIR}/generated/keywords) are not "
                      "those of the source; regenerate them with: ${GENERATOR} ${MANUAL_DIR}")
endif()

# Every page identical to its committed copy.
foreach(page keywords_table.rst compatibility_table.rst ${generated_entries})
  execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files ${WORK_DIR}/${page}
                          ${MANUAL_DIR}/${page} RESULT_VARIABLE differ)
  if (NOT differ EQUAL 0)
    message(FATAL_ERROR "${MANUAL_DIR}/${page} is out of date with the source; regenerate the "
                        "pages with: ${GENERATOR} ${MANUAL_DIR}")
  endif()
endforeach()
