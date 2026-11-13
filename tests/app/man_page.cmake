# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# The manual page of fepx (step 9.8, ROADMAP D101; CMake script, run with cmake -P by the test
# test_fepx_man_page): every option that "fepx --help" lists ("fepx --<option>") is described
# in the page (".B \-\-<option>"), and groff formats the page without a warning (when groff is
# found).
#
# Variables (-D): EXECUTABLE (fepx), MAN_PAGE (doc/man/fepx.1), GROFF (groff, or empty).

# The options of the help text.
execute_process(COMMAND ${EXECUTABLE} --help OUTPUT_VARIABLE help RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "fepx --help ended with status ${status}")
endif()
string(REGEX MATCHALL "fepx --[a-z-]+" options "${help}")
if(NOT options)
  message(FATAL_ERROR "fepx --help lists no option")
endif()

# Each one in the page, as an option entry.
file(READ ${MAN_PAGE} page)
foreach(option IN LISTS options)
  string(REPLACE "fepx --" "" name "${option}")
  string(FIND "${page}" ".B \\-\\-${name}\n" position)
  if(position EQUAL -1)
    message(FATAL_ERROR "${MAN_PAGE}: no entry for the option --${name} of fepx --help")
  endif()
endforeach()

# The page formats without a warning (-ww: all warnings; -z: no output).
if(GROFF)
  execute_process(COMMAND ${GROFF} -man -ww -z ${MAN_PAGE}
                  ERROR_VARIABLE warnings RESULT_VARIABLE status)
  if(NOT status EQUAL 0 OR NOT warnings STREQUAL "")
    message(FATAL_ERROR "groff: ${MAN_PAGE}: status ${status}\n${warnings}")
  endif()
endif()
