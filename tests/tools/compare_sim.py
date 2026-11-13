#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Compare two FEPX simulation result directories (.sim) field by field.

Source : FEPX regression harness (ROADMAP step 0.3). Development tool (CLAUDE.md §6: Python 3 +
         NumPy only for comparing outputs in tests). Tolerance policy: tests/regression/README.md.

Usage:
    compare_sim.py TEST_SIM REF_SIM [--rtol R] [--ftol F] [--atol A]
                   [--var-tol VAR=RTOL,FTOL ...] [--ignore GLOB ...] [--remove-mean VAR ...]
                   [--quiet]

TEST_SIM and REF_SIM are .sim directories (or any directories containing a `results/` tree).
Only files under `results/` are compared; `inputs/` and the `.sim` metadata file are not.

What is checked:
  1. The two `results/` trees contain the same set of files (apart from --ignore patterns).
     A missing or an extra file is a failure.
  2. Each pair of files has the same shape: same number of lines and of values per line.
     Comment lines (starting with '%', e.g. the headers of `forces` files) must be identical
     text (whitespace-insensitive).
  3. Every value agrees within the tolerance (all values on every line are compared):

         |test - ref| <= max(rtol * |ref|,  ftol * scale,  atol)

     where scale = max |ref| over the whole reference file (the field scale). The ftol term is
     the floor for components that are small compared with the rest of the field (e.g. shear
     stresses near zero in uniaxial loading), where a pure relative test is meaningless.
     A NaN or infinity in either file is a failure.

Variable-specific tolerances: --var-tol stress=1e-4,1e-6 applies to files of variable `stress`
(results/<entity>/stress/stress.stepN). The variable is the name of the directory holding the
file.

Mean translation: --remove-mean coo (repeatable) subtracts from the test values of files of that
variable, column by column, the mean over the lines of test - ref before comparing: a uniform
offset between the two runs is not a difference. Used for periodic runs (FEPX step 8.3b), whose
rigid translation FEPX fixes by the zero-mean fluctuation and v2 leaves to its solver
(dev/theory/addendum_periodic_mpc.md §5 V6): coordinates, displacements and velocities differ by
a uniform translation per step.

Exit status: 0 if all files agree, 1 if any check fails, 2 on a usage or read error.

Functions:
    list_result_files    relative paths of all files under <sim>/results
    read_values          parse one result file into a list of rows of floats
    compare_file         compare one pair of files, return a FileResult
    compare_sim          compare two .sim directories, return a list of FileResult
    main                 command-line interface
"""

import argparse
import fnmatch
import math
import os
import re
import sys
from dataclasses import dataclass, field

import numpy as np

# Default tolerances (rationale in tests/regression/README.md).
DEFAULT_RTOL = 1.0e-5
DEFAULT_FTOL = 1.0e-6
DEFAULT_ATOL = 0.0

# A Fortran Ew.d field drops the letter E when the exponent needs three digits: 1.0E-309 is
# written as ".1000000-308". Matches such tokens (mantissa, signed 3-digit exponent).
FORTRAN_E3 = re.compile(r"^([+-]?(?:\d+\.?\d*|\.\d+))([+-]\d{3})$")


@dataclass
class Tolerance:
    """Tolerances applied to one file.

    rtol  relative tolerance, per value (with respect to |ref|)
    ftol  tolerance relative to the field scale max|ref| of the file
    atol  absolute tolerance (in the units of the field)
    """

    rtol: float = DEFAULT_RTOL
    ftol: float = DEFAULT_FTOL
    atol: float = DEFAULT_ATOL


@dataclass
class FileResult:
    """Outcome of the comparison of one file.

    path        path relative to results/
    passed      True if the file agrees within the tolerance
    message     description of the failure, or a summary of the agreement
    max_abs     largest |test - ref| (NaN if not computed)
    max_scaled  largest |test - ref| / allowed difference (<= 1 means pass; NaN if not computed)
    """

    path: str
    passed: bool
    message: str
    max_abs: float = math.nan
    max_scaled: float = math.nan


@dataclass
class Options:
    """Comparison options.

    default   tolerance for all variables without an override
    per_var   tolerance overrides by variable name
    ignore    glob patterns (relative to results/) of files that are not compared
    remove_mean  variables whose mean test - ref offset is removed (--remove-mean)
    """

    default: Tolerance = field(default_factory=Tolerance)
    per_var: dict = field(default_factory=dict)
    ignore: list = field(default_factory=list)
    remove_mean: list = field(default_factory=list)


def list_result_files(sim_dir):
    """Return the sorted relative paths (POSIX separators) of all files under sim_dir/results."""
    results_dir = os.path.join(sim_dir, "results")
    if not os.path.isdir(results_dir):
        raise FileNotFoundError(f"no results directory: {results_dir}")
    paths = []
    for root, _dirs, files in os.walk(results_dir):
        for name in files:
            full = os.path.join(root, name)
            paths.append(os.path.relpath(full, results_dir).replace(os.sep, "/"))
    return sorted(paths)


def parse_number(token):
    """Convert one Fortran-formatted number to float (accepts D exponents and the E-less
    three-digit exponent form of Ew.d fields, e.g. '-.6982454-309')."""
    match = FORTRAN_E3.match(token)
    if match:
        return float(f"{match.group(1)}e{match.group(2)}")
    return float(token.replace("D", "E").replace("d", "e"))


def read_values(path):
    """Parse a result file into a list of rows.

    A numeric line becomes a list of floats; a comment line (first non-blank character '%')
    becomes its whitespace-normalized text (a str). Blank lines are kept as empty rows so
    that the line structure is compared too.
    """
    rows = []
    with open(path, "r", encoding="ascii") as handle:
        for line in handle:
            if line.lstrip().startswith("%"):
                rows.append(" ".join(line.split()))
            else:
                rows.append([parse_number(tok) for tok in line.split()])
    return rows


def variable_of(rel_path):
    """Return the variable name of a result file: the name of its parent directory."""
    parts = rel_path.split("/")
    return parts[-2] if len(parts) >= 2 else ""


def compare_file(rel_path, test_path, ref_path, tol, remove_mean=False):
    """Compare one pair of result files with the tolerance tol; return a FileResult.

    remove_mean: subtract the column means of test - ref from the test values first.
    """
    # Parse both files.
    try:
        test_rows = read_values(test_path)
        ref_rows = read_values(ref_path)
    except (OSError, ValueError, UnicodeDecodeError) as err:
        return FileResult(rel_path, False, f"cannot parse: {err}")

    # Shape: same number of lines, same number of values on each line.
    if len(test_rows) != len(ref_rows):
        return FileResult(rel_path, False,
                          f"line count differs: test {len(test_rows)}, ref {len(ref_rows)}")
    for i_line, (t_row, r_row) in enumerate(zip(test_rows, ref_rows), start=1):
        # Comment lines: identical text required.
        if isinstance(t_row, str) or isinstance(r_row, str):
            if t_row != r_row:
                return FileResult(rel_path, False, f"line {i_line}: comment/text differs: "
                                  f"test {t_row!r}, ref {r_row!r}")
            continue
        if len(t_row) != len(r_row):
            return FileResult(rel_path, False,
                              f"line {i_line}: value count differs: test {len(t_row)}, "
                              f"ref {len(r_row)}")

    # Flatten the numeric lines, keeping the (line, column) of each value for the report.
    numeric = [(i, t, r) for i, (t, r) in enumerate(zip(test_rows, ref_rows))
               if not isinstance(r, str)]
    # Mean translation removed (rows of equal length: node variables).
    if remove_mean and numeric and len({len(r) for _, _, r in numeric}) == 1:
        t_arr = np.array([t for _, t, _ in numeric], dtype=float)
        r_arr = np.array([r for _, _, r in numeric], dtype=float)
        t_arr = t_arr - (t_arr - r_arr).mean(axis=0)
        numeric = [(i, list(t_arr[k]), r) for k, (i, _, r) in enumerate(numeric)]
    test_vals = np.array([v for _, t, _ in numeric for v in t], dtype=float)
    ref_vals = np.array([v for _, _, r in numeric for v in r], dtype=float)
    if test_vals.size == 0:
        return FileResult(rel_path, True, "no numeric values in both", 0.0, 0.0)
    positions = [(i + 1, j + 1) for i, _, r in numeric for j in range(len(r))]

    # Non-finite values always fail.
    bad = ~(np.isfinite(test_vals) & np.isfinite(ref_vals))
    if bad.any():
        k = int(np.argmax(bad))
        line, col = positions[k]
        return FileResult(rel_path, False,
                          f"non-finite value at line {line}, column {col}: "
                          f"test {test_vals[k]!r}, ref {ref_vals[k]!r}")

    # Allowed difference per value: max(rtol |ref|, ftol scale, atol).
    scale = float(np.max(np.abs(ref_vals)))
    allowed = np.maximum(tol.rtol * np.abs(ref_vals), max(tol.ftol * scale, tol.atol))
    diff = np.abs(test_vals - ref_vals)
    max_abs = float(np.max(diff))

    # Ratio diff / allowed; where allowed is zero (exact zero field), any difference fails.
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(allowed > 0.0, diff / np.where(allowed > 0.0, allowed, 1.0),
                         np.where(diff > 0.0, np.inf, 0.0))
    k = int(np.argmax(ratio))
    max_scaled = float(ratio[k])
    if max_scaled > 1.0:
        line, col = positions[k]
        n_fail = int(np.count_nonzero(ratio > 1.0))
        return FileResult(rel_path, False,
                          f"{n_fail} of {ratio.size} values differ; worst at line {line}, "
                          f"column {col}: test {test_vals[k]:.9g}, ref {ref_vals[k]:.9g}, "
                          f"|diff| {diff[k]:.3g} > allowed {allowed[k]:.3g} "
                          f"(rtol {tol.rtol:g}, ftol {tol.ftol:g} x scale {scale:.3g}, "
                          f"atol {tol.atol:g})", max_abs, max_scaled)
    return FileResult(rel_path, True, f"max |diff| {max_abs:.3g}, max diff/allowed "
                      f"{max_scaled:.3g}", max_abs, max_scaled)


def is_ignored(rel_path, patterns):
    """Return True if rel_path matches one of the glob patterns."""
    return any(fnmatch.fnmatch(rel_path, pat) for pat in patterns)


def compare_sim(test_sim, ref_sim, options):
    """Compare all result files of two .sim directories; return a list of FileResult."""
    # File sets, without ignored files.
    test_files = [p for p in list_result_files(test_sim) if not is_ignored(p, options.ignore)]
    ref_files = [p for p in list_result_files(ref_sim) if not is_ignored(p, options.ignore)]
    test_set, ref_set = set(test_files), set(ref_files)

    results = []
    # Missing and extra files.
    for path in sorted(ref_set - test_set):
        results.append(FileResult(path, False, "missing in test"))
    for path in sorted(test_set - ref_set):
        results.append(FileResult(path, False, "extra in test (not in reference)"))

    # Common files, with the tolerance of their variable.
    for path in sorted(test_set & ref_set):
        tol = options.per_var.get(variable_of(path), options.default)
        results.append(compare_file(path, os.path.join(test_sim, "results", path),
                                    os.path.join(ref_sim, "results", path), tol,
                                    variable_of(path) in options.remove_mean))
    return results


def parse_var_tol(text):
    """Parse VAR=RTOL,FTOL[,ATOL] into (var, Tolerance)."""
    try:
        var, values = text.split("=", 1)
        numbers = [float(x) for x in values.split(",")]
        if not var or len(numbers) not in (2, 3):
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected VAR=RTOL,FTOL[,ATOL], got {text!r}")
    return var, Tolerance(*numbers)


def main(argv=None):
    """Command-line interface; returns the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("test_sim", help="simulation.sim directory to check")
    parser.add_argument("ref_sim", help="reference .sim directory")
    parser.add_argument("--rtol", type=float, default=DEFAULT_RTOL, help="relative tolerance")
    parser.add_argument("--ftol", type=float, default=DEFAULT_FTOL,
                        help="tolerance relative to the field scale max|ref| of each file")
    parser.add_argument("--atol", type=float, default=DEFAULT_ATOL, help="absolute tolerance")
    parser.add_argument("--var-tol", type=parse_var_tol, action="append", default=[],
                        metavar="VAR=RTOL,FTOL[,ATOL]", help="tolerance for one variable")
    parser.add_argument("--ignore", action="append", default=[], metavar="GLOB",
                        help="glob (relative to results/) of files not to compare")
    parser.add_argument("--remove-mean", action="append", default=[], metavar="VAR",
                        help="variable whose mean test - ref offset is removed (periodic runs)")
    parser.add_argument("--quiet", action="store_true", help="print only failures and summary")
    args = parser.parse_args(argv)

    options = Options(Tolerance(args.rtol, args.ftol, args.atol), dict(args.var_tol), args.ignore,
                      args.remove_mean)
    try:
        results = compare_sim(args.test_sim, args.ref_sim, options)
    except FileNotFoundError as err:
        print(f"compare_sim: {err}", file=sys.stderr)
        return 2

    # Report.
    n_fail = 0
    for res in results:
        if not res.passed:
            n_fail += 1
            print(f"FAIL  {res.path}: {res.message}")
        elif not args.quiet:
            print(f"ok    {res.path}: {res.message}")
    worst = max((r.max_scaled for r in results if not math.isnan(r.max_scaled)), default=0.0)
    print(f"compare_sim: {len(results)} files, {n_fail} failed, "
          f"worst diff/allowed {worst:.3g}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
