# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Reading FEPX's terminal output in the tests (step 9.7, ROADMAP D94).

The log's messages are wrapped at 72 columns (src/core/log_mod.f90): a longer message continues
on lines that start with its prefix ("Info   :", "Warning:", "Error  :") followed by two or more
blanks and no depth marker ("- ", "> ", ". ") or file marker ("[i] ", "[o] "). joined_lines
returns the log's lines with every continuation joined to its message with one blank (the rule
of tests/core/unwrap_log.cmake). Development tool (CLAUDE.md §6).
"""

PREFIXES = ("Info   :", "Warning:", "Error  :")
MARKERS = ("- ", "> ", ". ", "[i] ", "[o] ")


def joined_lines(text):
    """The lines of a log, each wrapped message joined into one line."""
    lines = []
    previous = None
    for line in text.splitlines():
        prefix = line[:8] if line[:8] in PREFIXES else None
        if prefix is not None and prefix == previous and line[9:11] == "  ":
            rest = line[9:].lstrip(" ")
            if rest and not rest.startswith(MARKERS):
                lines[-1] += " " + rest
                continue
        lines.append(line)
        previous = prefix
    return lines
