#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check that every option-compatibility rule is used and tested (step 9.9, ROADMAP D104).

Usage: check_rules_used.py <rules file> <source directory>

The rules file is written by fepx_keywords --rules: one line "<name>|<message>" per rule of
src/io/cfg_rules_mod.f90, the message being that of the subject "subject". For every rule:
  - the configuration builder (src/io/config_mod.f90) names its index rule_<name>;
  - an error test expects its message: the relation and the condition ("used only with ...",
    "cannot be used with ...", "ignored with ...", up to the parenthesis of the reason) appear
    in an expected message of a test definition (tests/**/CMakeLists.txt, tests/**/*.cmake),
    where the tests write "(" and ")" as "[(]" and "[)]" and may write any of . ; : " ' as ".".
Prints every rule that fails and exits with status 1 if any does.

Source: FEPX (ROADMAP step 9.9, D104). A development tool of the test suite (CLAUDE.md §6).
"""

import glob
import os
import re
import sys


def normalized(text):
    """The text with the tests' escapes undone and the characters they may write as "." made ".".

    text  (str) a message or a test definition
    """
    text = text.replace("[(]", "(").replace("[)]", ")")
    return re.sub(r"[.;:\"']", ".", text)


def main():
    """Check the rules; return the exit status."""
    rules_path, source_dir = sys.argv[1], sys.argv[2]
    with open(rules_path, encoding="utf-8") as handle:
        rules = [line.rstrip("\n").split("|", 1) for line in handle if line.strip()]

    # The builder's source and the test definitions.
    with open(os.path.join(source_dir, "src", "io", "config_mod.f90"), encoding="utf-8") as handle:
        builder = handle.read()
    test_files = glob.glob(os.path.join(source_dir, "tests", "**", "CMakeLists.txt"),
                           recursive=True)
    test_files += glob.glob(os.path.join(source_dir, "tests", "**", "*.cmake"), recursive=True)
    tests = ""
    for path in sorted(test_files):
        with open(path, encoding="utf-8") as handle:
            tests += normalized(handle.read())

    # Every rule: used by its index, its message expected.
    failures = []
    for name, message in rules:
        if not re.search(r"\brule_" + re.escape(name) + r"\b", builder):
            failures.append(f"rule {name}: rule_{name} is not used in src/io/config_mod.f90")
        relation = re.sub(r"^subject (is |are )?", "", message).split(" (", 1)[0]
        if normalized(relation) not in tests:
            failures.append(f"rule {name}: no test expects the message '... {relation}'")
    for failure in failures:
        print(failure)
    print(f"{len(rules)} rules checked, {len(failures)} failures")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
