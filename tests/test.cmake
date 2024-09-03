# This file is part of the FEPX software package.
# Copyright (C) 1996-2021, DPLab, ACME Lab.
# See the COPYING file in the top-level directory.

if (DEFINED post_command)
  execute_process(COMMAND ${post_command})
endif()

if (DEFINED pre_command)
  execute_process(COMMAND ${pre_command})
endif()

execute_process(COMMAND mpirun -np 2 ${test_prog} RESULT_VARIABLE RESVAR)

if(RESVAR)
  message(FATAL_ERROR "Test failed")
endif()

if ("${test_mode}" MATCHES "Normal" AND NOT "${test_mode_force_minimal}" EQUAL 1)

  file(GLOB bak_files *.bak *~ FAILED)
  foreach(bak_file ${bak_files})
    file (REMOVE ${bak_file})
  endforeach()

  if (DEFINED ref_dirs)
    file(GLOB_RECURSE ref_files ref.sim/**/*)
  endif()


  foreach(ref_file ${ref_files})
    string(REPLACE "ref.sim" "simulation.sim" test_file ${ref_file})
    if (NOT "${test_file}" MATCHES ".png")
      if ("${test_diff}" MATCHES "Hard")
        execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files
                        ${test_file} ${ref_file}
                        RESULT_VARIABLE RESVAR)
      else()
        execute_process(COMMAND ${test_prog} --diff
                        ${test_file} ${ref_file}
                        RESULT_VARIABLE RESVAR)
      endif()
    else()
      execute_process(COMMAND compare -metric AE ${test_file} ${ref_file} NULL:
                      RESULT_VARIABLE RESVAR)
    endif()

    if(RESVAR)
      file(RENAME ${test_file} ${test_file}.bak)
      message(FATAL_ERROR "Test failed - files differ")
      file (REMOVE ${test_file})
    endif()
  endforeach()

elseif ("${test_mode}" MATCHES "Writing")
  file(REMOVE_RECURSE ref.sim)
  file(RENAME simulation.sim ref.sim)
endif()
