#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Compare FEPX's terminal output with its reference after masking the variable fields (FEPX
step 9.7, ROADMAP D94; dev/plan/terminal_output.md §3.8, the golden-file test).

Source : FEPX. Development tool (CLAUDE.md §6: Python 3 for comparing outputs in tests).

The whole standard output of a run is compared line by line with the reference, so that the
structure (banner, setup summary, progress lines, closing summary), the wording, the depths and
the line breaks are checked. Masked first, in the output and the reference alike, are the fields
that differ between builds, machines and runs:
  - the version line (git hash and date) and the build line (compiler, MPI, HDF5);
  - every run of global-iteration lines ("> Iteration <k> (<SA|NR>): R = ... (<n> CG
    iterations)") becomes one line "> Iteration <masked>": the number of iterations, their
    methods, R and the CG counts depend on the floating-point details of the build (optimization
    level, MPI library; results agree to solver tolerance only);
  - "Converged in <n> iterations (<t> s)", the closing summary's global-iteration and CG counts;
  - every wall time ("<t> s" in the closing summary's table, in the columns of the time, so
    that the times' alignment stays checked; "(<t> s)").
Everything else must be identical: the counts of increments, the times and lengths of the
increments, the input summary, the partition statistics, the warnings issued.

Usage:
  check_terminal_output.py OUTPUT REFERENCE
With the environment variable FEPX_REGRESSION_UPDATE set (to 1 or all), the masked output is
written to REFERENCE instead (tests/regression/README.md: only with a documented reason).
Exit status 0 when the masked output equals the reference, 1 otherwise (with a unified diff).
"""

import difflib
import os
import re
import sys


ITERATION = re.compile(r"^(Info   :\s+> )Iteration \d+ \((SA|NR)\): R = ")


def masked(text):
    """The output's lines with the variable fields masked (see the module docstring)."""
    lines = []
    for line in text.splitlines():
        if line.startswith("Info   : Version "):
            line = "Info   : Version <masked>"
        elif line.startswith("Info   : Built with: "):
            line = "Info   : Built with: <masked>"
        match = ITERATION.match(line)
        if match:
            line = match.group(1) + "Iteration <masked>"
            if lines and lines[-1] == line:
                continue
        line = re.sub(r"Converged in \d+ iterations? \(\d+\.\d+ s\)",
                      "Converged in <n> iterations (<t> s)", line)
        line = re.sub(r"Global iterations: \d+ \(\d+ NR\)", "Global iterations: <n> (<n> NR)",
                      line)
        line = re.sub(r"CG iterations: \d+$", "CG iterations: <n>", line)
        # (A table time masked in its own columns: the alignment of the "s" stays checked.)
        line = re.sub(r":(\s+)(\d+\.\d{3}) s$",
                      lambda m: ":" + " " * (len(m.group(1)) + len(m.group(2)) - 3) + "<t> s",
                      line)
        lines.append(line)
    return lines


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    output_path, reference_path = sys.argv[1:]
    with open(output_path, encoding="utf-8", errors="replace") as stream:
        output = masked(stream.read())
    if os.environ.get("FEPX_REGRESSION_UPDATE", "") in ("1", "all"):
        with open(reference_path, "w", encoding="utf-8") as stream:
            stream.write("\n".join(output) + "\n")
        print(f"reference written: {reference_path}")
        return 0
    with open(reference_path, encoding="utf-8") as stream:
        reference = stream.read().splitlines()
    if output == reference:
        print(f"terminal output equals the reference ({len(output)} lines, masked)")
        return 0
    print("FAIL  the terminal output differs from the reference (masked):")
    sys.stdout.writelines(line + "\n" for line in difflib.unified_diff(
        reference, output, "reference", "output", lineterm=""))
    return 1


if __name__ == "__main__":
    sys.exit(main())
