#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check that every user-facing error and warning message has an entry in the manual.

Usage: check_messages.py <source directory> [--list]

Every call of fatal_error, cfg_error, line_error, text_error and log_warning in src/ (the programs' own
usage messages excepted: src/app/fepx_keywords.f90, src/app/fepx_slip_systems.f90) gives a
message; its *key* is the longest literal fragment of the message (the routine name, the first
argument of fatal_error, left out). When the call passes a message built in a variable, the key
is taken from the literals of the last "write (<variable>, ...)" before the call in the same
routine. The key, with runs of blanks made one blank, must appear in
doc/messages.rst. Not checked:
  - calls of internal_error (internal consistency checks, not caused by the input: ROADMAP
    D104);
  - messages of the option-compatibility rules (cfg_rule_message), documented by the
    compatibility table generated from src/io/cfg_rules_mod.f90;
  - messages whose key is shorter than 8 characters (a message given whole by the caller: its
    own call is checked).
With --list, prints every call and its key instead (for writing the manual's entries).
Prints every message without an entry and exits with status 1 if any.

Source: FEPX (ROADMAP step 9.9, D104). A development tool of the test suite (CLAUDE.md §6).
"""

import glob
import os
import re
import sys

CALL = re.compile(r"\bcall (fatal_error|cfg_error|line_error|text_error|log_warning)\s*\(")
# (internal_error calls are not matched: internal consistency checks.)
# A literal in single or double quotes (doubled quotes inside): the text of either group.
LITERAL_TOKEN = re.compile(r"'((?:[^']|'')*)'|\"((?:[^\"]|\"\")*)\"")
EXCLUDED = ("src/app/fepx_keywords.f90", "src/app/fepx_slip_systems.f90")
MIN_KEY = 8
# A format of a write statement, e.g. "(a,i0,a,es13.6)".
FORMAT = re.compile(r"^\((\s*(\d*[aif]\d*|es\d+\.\d+|f\d+\.\d+|g0|i0|a)\s*,?)+\)$")


def literals_of(text):
    """The literals of a statement, in order, each with its doubled quotes made single."""
    result = []
    for match in LITERAL_TOKEN.finditer(text):
        if match.group(1) is not None:
            result.append(match.group(1).replace("''", "'"))
        else:
            result.append(match.group(2).replace('""', '"'))
    return result


def statements(path):
    """The statements of a Fortran source file, continuation lines joined, comments removed.

    path        (str) the source file
    statements  (list of (int, str)) first line number and text of each statement
    """
    result = []
    text = ""
    first = 0
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            # Remove a comment that is not inside a literal: scan the quotes.
            code, quote = "", ""
            for char in line.rstrip("\n"):
                if char in "'\"" and quote in ("", char):
                    quote = "" if quote == char else char
                if char == "!" and not quote:
                    break
                code += char
            code = code.strip()
            if not text:
                first = number
            if code.startswith("&"):
                code = code[1:]
            continued = code.endswith("&")
            text += " " + (code[:-1] if continued else code)
            if not continued:
                if text.strip():
                    result.append((first, text.strip()))
                text = ""
    return result


def key_of(literals):
    """The longest literal fragment (blanks collapsed), formats of
    write statements left out."""
    fragments = [re.sub(r"\s+", " ", item).strip() for item in literals
                 if not FORMAT.match(item)]
    return max(fragments + [""], key=len)


def messages(source_dir):
    """Every checked call: (file, line, key) (key "" when none could be found)."""
    found = []
    for path in sorted(glob.glob(os.path.join(source_dir, "src", "**", "*.f90"),
                                 recursive=True)):
        relative = os.path.relpath(path, source_dir)
        if relative in EXCLUDED:
            continue
        last_write = {}
        for line, text in statements(path):
            # A routine boundary forgets the messages built so far.
            if re.match(r"(end )?(subroutine|function)\b|.*\b(subroutine|function) \w+\s*\(",
                        text):
                last_write = {}
            write = re.match(r"write\s*\(\s*(\w+)\s*,", text)
            if write:
                last_write[write.group(1)] = literals_of(text)
            call = CALL.search(text)
            if not call or "cfg_rule_message" in text:
                continue
            arguments = text[call.end():]
            literals = literals_of(arguments)
            if call.group(1) == "fatal_error" and literals and arguments.lstrip().startswith("'"):
                literals = literals[1:]
            key = key_of(literals)
            if len(key) < MIN_KEY:
                # A message built in a variable: the literals of its write.
                for name, written in last_write.items():
                    if re.search(r"\b" + name + r"\b", arguments):
                        key = key_of(written)
                        break
            found.append((relative, line, key))
    return found


def main():
    """Check (or list) the messages; return the exit status."""
    source_dir = sys.argv[1]
    found = messages(source_dir)
    if "--list" in sys.argv[2:]:
        for relative, line, key in found:
            print(f"{relative}:{line}: {key}")
        return 0
    with open(os.path.join(source_dir, "doc", "messages.rst"),
              encoding="utf-8") as handle:
        manual = re.sub(r"\s+", " ", handle.read().replace("``", ""))
    missing = [(relative, line, key) for relative, line, key in found
               if len(key) >= MIN_KEY and key not in manual]
    for relative, line, key in missing:
        print(f"{relative}:{line}: no entry in doc/messages.rst for '{key}'")
    print(f"{len(found)} messages checked, {len(missing)} without an entry")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
