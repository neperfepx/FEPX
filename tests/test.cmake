# This file is part of the FEPX software package.
# Copyright (C) 1996-2021, DPLab, ACME Lab.
# See the COPYING file in the top-level directory.

execute_process(COMMAND mpirun -np 2 ${test_prog} RESULT_VARIABLE RESVAR)

if(RESVAR)
  message(FATAL_ERROR "Test failed")
endif()

if ("${test_mode}" MATCHES "Normal" AND NOT "${test_mode_force_minimal}" EQUAL 1)
  file(REMOVE FAILED)
  file(GLOB ref_files refpost.*)
  foreach(ref_file ${ref_files})
    string(REPLACE "refpost." "post." test_file ${ref_file})
    execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files
                    ${test_file} ${ref_file}
                    RESULT_VARIABLE RESVAR)
    if(RESVAR)
      file(WRITE FAILED "Test failed - files ${ref_file} and ${test_file} differ")
      message(FATAL_ERROR "Test failed - files ${ref_file} and ${test_file} differ")
    else()
      file (REMOVE ${test_file})
    endif()
  endforeach()

elseif ("${test_mode}" MATCHES "Writing")
  file(GLOB test_files post.*)
  foreach(test_file ${test_files})
    string(REPLACE "post." "refpost." ref_file ${test_file})
    file(RENAME ${test_file} ${ref_file})
  endforeach()
endif()
