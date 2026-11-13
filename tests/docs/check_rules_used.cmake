# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Check that every option-compatibility rule (src/io/cfg_rules_mod.f90, step 9.9, ROADMAP D104)
# is used by the configuration builder and that its message is expected by an error test: the
# rules are listed by fepx_keywords --rules, then checked by check_rules_used.py.
#
# Variables: GENERATOR (fepx_keywords executable), PYTHON (Python 3), SCRIPT
# (check_rules_used.py), SOURCE_DIR (the source tree), RULES_FILE (scratch file in the build
# tree).

# The list of the rules.
execute_process(COMMAND ${GENERATOR} --rules ${RULES_FILE} RESULT_VARIABLE status)
if (NOT status EQUAL 0)
  message(FATAL_ERROR "${GENERATOR} --rules failed (status ${status})")
endif()

# The check.
execute_process(COMMAND ${PYTHON} ${SCRIPT} ${RULES_FILE} ${SOURCE_DIR} RESULT_VARIABLE status)
if (NOT status EQUAL 0)
  message(FATAL_ERROR "${SCRIPT} failed (status ${status})")
endif()
