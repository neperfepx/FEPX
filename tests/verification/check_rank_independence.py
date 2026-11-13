# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check that v2-fixed runs of one case on different rank counts give the same element stresses.

Usage: check_rank_independence.py --rtol=<tol> <run dir 1> <run dir 2> [...]

Every printed step of results/elts/stress in each run directory's simulation.sim is compared
with the first run; the largest difference relative to the largest |stress| of the step must
not exceed rtol. Prints one line per run. Development tool (CLAUDE.md §6); CHANGES.md C36.
"""
import argparse
import glob
import os
import sys

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--rtol', type=float, required=True)
    parser.add_argument('runs', nargs='+')
    args = parser.parse_args()

    ref = args.runs[0]
    pattern = os.path.join(ref, 'simulation.sim', 'results', 'elts', 'stress', 'stress.step*')
    ref_files = sorted(glob.glob(pattern))
    if not ref_files:
        print(f'no stress files in {ref}')
        return 1
    failed = False
    for run in args.runs[1:]:
        worst = 0.0
        for ref_file in ref_files:
            other = os.path.join(run, os.path.relpath(ref_file, ref))
            a = np.loadtxt(ref_file)
            b = np.loadtxt(other)
            scale = np.abs(a).max()
            if scale > 0.0:
                worst = max(worst, np.abs(a - b).max() / scale)
        ok = worst <= args.rtol
        failed = failed or not ok
        print(f"{'ok  ' if ok else 'FAIL'} {os.path.basename(run)} vs {os.path.basename(ref)}: "
              f'largest relative stress difference {worst:.3e} (tolerance {args.rtol:.1e})')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
