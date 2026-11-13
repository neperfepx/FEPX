#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check that every configuration input shown in the manual is a file the test suite runs.

Usage: check_manual_inputs.py <manual directory>

A literal block of the manual (after "::", or a code-block or code directive) is an *input*
when at least two of its lines start with a keyword of FEPX (the keywords of the generated
entries, <manual>/generated/keywords/<keyword>.inc) and none of its lines has a mark of a
syntax template ("<", "[ " or "..."). Inputs must be included from tested files
(literalinclude, which this check does not read: ROADMAP D104), not written in the page; an
input that is deliberately not run (an FEPX v2 input, a variant for the reader) is preceded,
within the five lines before it, by a comment ".. Not run: <reason>". Prints every other input
block and exits with status 1 if any.

Source: FEPX (ROADMAP step 9.9, D104). A development tool of the test suite (CLAUDE.md §6).
"""

import glob
import os
import re
import sys

GENERATED = ("keywords_table.rst", "compatibility_table.rst", "slip_systems_table.rst")


def literal_blocks(lines):
    """The literal blocks of a page: (line number of the introducing line, body lines)."""
    blocks = []
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        introduces = (stripped.endswith("::") and not stripped.startswith("..")) or \
            re.match(r"\.\. (code-block|code)::", stripped)
        if introduces:
            indent = len(line) - len(line.lstrip())
            j = i + 1
            body = []
            while j < len(lines) and (not lines[j].strip()
                                      or len(lines[j]) - len(lines[j].lstrip()) > indent):
                # Directive options (":language: none") are not part of the body.
                if lines[j].strip() and not lines[j].strip().startswith(":"):
                    body.append(lines[j].strip())
                j += 1
            blocks.append((i + 1, body))
            i = j
        else:
            i += 1
    return blocks


def main():
    """Check the pages; return the exit status."""
    manual_dir = sys.argv[1]
    keywords = {os.path.basename(path)[:-len(".inc")] for path in
                glob.glob(os.path.join(manual_dir, "generated", "keywords", "*.inc"))}
    failures = 0
    for path in sorted(glob.glob(os.path.join(manual_dir, "**", "*.rst"), recursive=True)):
        if os.path.basename(path) in GENERATED or os.sep + "generated" + os.sep in path:
            continue
        with open(path, encoding="utf-8") as handle:
            lines = handle.read().split("\n")
        for number, body in literal_blocks(lines):
            keyword_lines = sum(1 for line in body if line.split()[0] in keywords)
            template = any("<" in line or "[ " in line or "..." in line for line in body)
            if keyword_lines < 2 or template:
                continue
            before = lines[max(0, number - 6):number]
            if any(line.strip().startswith(".. Not run:") for line in before):
                continue
            print(f"{os.path.relpath(path, manual_dir)}:{number}: a configuration input written "
                  "in the page: include it from a tested file (literalinclude)")
            failures += 1
    print(f"{failures} input(s) written in the pages")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
