# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Source archives of FEPX (step 9.8, ROADMAP D101; CMake script, run with cmake -P by the targets
# dist and dist_testdata of src/cmake/Packaging.cmake).
#
# Variables (-D):
#   SOURCE_DIR  the source tree (a git repository)
#   OUTPUT_DIR  where the archive is written
#   NAME        the archive's name and top directory, e.g. fepx-3.0.0-dev
#   PART        source: <NAME>.tar.gz, the tree without the regression reference data
#               (tests/regression/reference/) and without the files marked export-ignore in
#               .gitattributes; testdata: <NAME>-testdata.tar.gz, the reference data alone, with
#               the same top directory, so that it unpacks into the source tree
#
# The archive is that of HEAD, the committed state (git archive): uncommitted changes are not in
# it, and a warning says so. src/cmake/git_archive_info.txt gets the commit's hash and date
# (export-subst), so that a build from the archive shows them (cmake/BuildInfo.cmake).

find_program(GIT_EXECUTABLE git)
if(NOT GIT_EXECUTABLE)
  message(FATAL_ERROR "dist: git is needed to make the source archives.")
endif()

# Uncommitted changes of tracked files are not archived.
execute_process(COMMAND ${GIT_EXECUTABLE} -C ${SOURCE_DIR} status --porcelain --untracked-files=no
                OUTPUT_VARIABLE changes RESULT_VARIABLE status OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "dist: ${SOURCE_DIR} is not a git repository.")
endif()
if(NOT changes STREQUAL "")
  message(WARNING "dist: the source tree has uncommitted changes; the archive is that of HEAD.")
endif()

# The archive and the paths it holds (git pathspecs).
set(references tests/regression/reference)
if(PART STREQUAL "source")
  set(archive ${OUTPUT_DIR}/${NAME}.tar.gz)
  set(paths . ":(exclude)${references}")
elseif(PART STREQUAL "testdata")
  set(archive ${OUTPUT_DIR}/${NAME}-testdata.tar.gz)
  set(paths ${references})
else()
  message(FATAL_ERROR "dist: PART '${PART}' is not source or testdata.")
endif()
execute_process(COMMAND ${GIT_EXECUTABLE} -C ${SOURCE_DIR} archive --format=tar.gz
                        --prefix=${NAME}/ --output=${archive} HEAD -- ${paths}
                RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "dist: git archive failed.")
endif()
message(STATUS "dist: ${archive}")
