# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Installation and source archives (step 9.8, ROADMAP D101; no Debian package since D122).
#
# Installed (cmake --install, GNUInstallDirs layout under CMAKE_INSTALL_PREFIX):
#   bin/fepx                              the executable (src/CMakeLists.txt)
#   share/man/man1/fepx.1[.gz]            the manual page (doc/man/fepx.1; gzipped when gzip
#                                         is found)
#   share/doc/fepx/README.md, COPYING     the README and the licence
#   share/doc/fepx/examples/<case>/       the manual's tutorial cases (tutorials/
#                                         tutorials.txt): simulation.cfg (.phase, .ori) and the
#                                         mesh as simulation.msh, ready to run; and the
#                                         post-processing tutorial's script
# Not installed: fepx_keywords and fepx_slip_systems (the manual's tables), the libraries,
# the HTML manual (online).
#
# The conda packages are built from this installation (conda/recipe, dev/packaging.md, D121).
#
# Target userinstall (as Neper's make userinstall): builds fepx and installs it as above under
# $HOME/.local (whatever CMAKE_INSTALL_PREFIX is), without administrator rights
# (src/cmake/UserInstall.cmake).
#
# Source archives (git archive of HEAD, the committed state; src/cmake/Dist.cmake):
#   target dist           fepx-<version>.tar.gz, without the regression reference data and
#                         without the files marked export-ignore in .gitattributes
#   target dist_testdata  fepx-<version>-testdata.tar.gz, the reference data alone
#                         (tests/regression/reference/), unpacked over the source tree
#
# Variables used: FEPX_PACKAGE_NAME (CMakeLists.txt), FEPX_VERSION_SUFFIX.

# The manual page, gzipped when gzip is found (-n: no name or time stamp, so
# that the same page gives the same file).
set(fepx_man_page ${FEPX_ROOT_DIR}/doc/man/fepx.1)
find_program(FEPX_GZIP_PROGRAM gzip)
if(FEPX_GZIP_PROGRAM)
  add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/man/fepx.1.gz
                     COMMAND ${CMAKE_COMMAND} -E copy ${fepx_man_page} ${PROJECT_BINARY_DIR}/man/fepx.1
                     COMMAND ${FEPX_GZIP_PROGRAM} -n -9 -f ${PROJECT_BINARY_DIR}/man/fepx.1
                     DEPENDS ${fepx_man_page}
                     COMMENT "Manual page fepx.1.gz" VERBATIM)
  add_custom_target(fepx_man_page ALL DEPENDS ${PROJECT_BINARY_DIR}/man/fepx.1.gz)
  install(FILES ${PROJECT_BINARY_DIR}/man/fepx.1.gz DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
else()
  install(FILES ${fepx_man_page} DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
endif()

# The documentation files.
install(FILES ${FEPX_ROOT_DIR}/README.md ${FEPX_ROOT_DIR}/COPYING
        DESTINATION ${CMAKE_INSTALL_DOCDIR})

# The tutorial cases, one directory each, from the list of the manual's tutorials (its columns:
# case, inputs directory, mesh, orientation file, the case it continues, ranks, tier).
set(fepx_tutorials_dir ${FEPX_ROOT_DIR}/tutorials)
set(fepx_meshes_dir ${FEPX_ROOT_DIR}/tests/regression/cases/meshes)
set(fepx_examples_dir ${CMAKE_INSTALL_DOCDIR}/examples)
file(STRINGS ${fepx_tutorials_dir}/tutorials.txt fepx_tutorial_lines)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${fepx_tutorials_dir}/tutorials.txt)
foreach(line IN LISTS fepx_tutorial_lines)
  # Comments and empty lines skipped; the columns as a list.
  string(REGEX REPLACE "#.*$" "" line "${line}")
  string(STRIP "${line}" line)
  if(line STREQUAL "")
    continue()
  endif()
  string(REGEX REPLACE "[ \t]+" ";" columns "${line}")
  list(GET columns 0 case_name)
  list(GET columns 1 inputs)
  list(GET columns 2 mesh)
  list(GET columns 3 ori)
  # The configuration and the phase file of the inputs directory (not its other orientation
  # files, which belong to other cases).
  file(GLOB case_files ${fepx_tutorials_dir}/${inputs}/simulation.cfg
                       ${fepx_tutorials_dir}/${inputs}/simulation.phase)
  install(FILES ${case_files} DESTINATION ${fepx_examples_dir}/${case_name})
  install(FILES ${fepx_meshes_dir}/${mesh} DESTINATION ${fepx_examples_dir}/${case_name}
          RENAME simulation.msh)
  # The orientation file: in the inputs directory ("inputs:<file>") or with the meshes.
  if(ori MATCHES "^inputs:(.+)$")
    install(FILES ${fepx_tutorials_dir}/${inputs}/${CMAKE_MATCH_1}
            DESTINATION ${fepx_examples_dir}/${case_name} RENAME simulation.ori)
  elseif(NOT ori STREQUAL "-")
    install(FILES ${fepx_meshes_dir}/${ori} DESTINATION ${fepx_examples_dir}/${case_name}
            RENAME simulation.ori)
  endif()
endforeach()
install(FILES ${fepx_tutorials_dir}/postprocessing/postprocess.py
        DESTINATION ${fepx_examples_dir}/postprocessing)

# The installation for the user, under $HOME/.local (cmake/UserInstall.cmake; HOME read when the
# target runs). cmake --install does not build: the target depends on what is installed and built
# (the executable and, when gzipped, the manual page).
add_custom_target(userinstall
                  COMMAND ${CMAKE_COMMAND} -DBINARY_DIR=${PROJECT_BINARY_DIR}
                          -P ${PROJECT_SOURCE_DIR}/cmake/UserInstall.cmake
                  COMMENT "Installing FEPX into $HOME/.local" USES_TERMINAL VERBATIM)
add_dependencies(userinstall fepx)
if(TARGET fepx_man_page)
  add_dependencies(userinstall fepx_man_page)
endif()

# The source archives (cmake/Dist.cmake; git archive of HEAD).
set(fepx_dist_command ${CMAKE_COMMAND} -DSOURCE_DIR=${FEPX_ROOT_DIR}
    -DOUTPUT_DIR=${PROJECT_BINARY_DIR}
    -DNAME=${FEPX_PACKAGE_NAME}-${PROJECT_VERSION}${FEPX_VERSION_SUFFIX})
add_custom_target(dist COMMAND ${fepx_dist_command} -DPART=source
                       -P ${PROJECT_SOURCE_DIR}/cmake/Dist.cmake
                  COMMENT "Source tarball (git archive of HEAD)" VERBATIM)
add_custom_target(dist_testdata COMMAND ${fepx_dist_command} -DPART=testdata
                                -P ${PROJECT_SOURCE_DIR}/cmake/Dist.cmake
                  COMMENT "Test-data tarball (git archive of HEAD)" VERBATIM)
