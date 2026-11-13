# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# The joining of the log's wrapped messages (step 9.7; included by check_run.cmake and the other
# test scripts that match FEPX's terminal output).
#
# Join the continuation lines of the log's messages: a line "<prefix><two or more blanks><text>"
# whose text does not start with a depth marker ("- ", "> ", ". ") or a file marker ("[i] ",
# "[o] "), after a line with the same prefix ("Info   :", "Warning:", "Error  :"), is appended
# to that line with one blank.
function(fepx_unwrap_log text out_var)
  # Protect the semicolons (list separators), then split into lines.
  string(REPLACE ";" "<fepx-semicolon>" text "${text}")
  string(REPLACE "\n" ";" lines "${text}")
  set(result "")
  set(previous_prefix "")
  set(first TRUE)
  foreach(line IN LISTS lines)
    set(joined FALSE)
    if(line MATCHES "^(Info   :|Warning:|Error  :) *(.*)$")
      set(prefix "${CMAKE_MATCH_1}")
      set(rest "${CMAKE_MATCH_2}")
      if(line MATCHES "^(Info   :|Warning:|Error  :)   " AND prefix STREQUAL previous_prefix AND
         NOT rest MATCHES "^([->.] |\\[[io]\\] )")
        string(APPEND result " ${rest}")
        set(joined TRUE)
      endif()
      set(previous_prefix "${prefix}")
    else()
      set(previous_prefix "")
    endif()
    if(NOT joined)
      if(first)
        set(result "${line}")
        set(first FALSE)
      else()
        string(APPEND result "\n${line}")
      endif()
    endif()
  endforeach()
  string(REPLACE "<fepx-semicolon>" ";" result "${result}")
  set(${out_var} "${result}" PARENT_SCOPE)
endfunction()

# The log made comparable: its messages joined (fepx_unwrap_log), every run of blanks one blank.
function(fepx_normalize_log text out_var)
  fepx_unwrap_log("${text}" joined)
  string(REGEX REPLACE "  +" " " joined "${joined}")
  set(${out_var} "${joined}" PARENT_SCOPE)
endfunction()
