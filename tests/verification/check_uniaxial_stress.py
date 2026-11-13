#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check that the stress of a homogeneous uniaxial-stress run is uniaxial in every element.

Source : FEPX verification tests (ROADMAP step 0.4, CHANGES.md C11). Development tool
         (CLAUDE.md §6).

A homogeneous single crystal loaded along z with traction-free lateral faces (uniaxial_minimal)
must carry uniaxial stress: in every element, every stress component other than sigma_zz must be
below --tol * |sigma_zz| at step --step (which must exist).

Usage: check_uniaxial_stress.py SIM_DIR --step N [--tol T]
Exit status 0 if the check passes, 1 otherwise.
"""

import argparse
import os
import sys

import numpy as np


def main(argv=None):
    """Command-line interface; returns the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("sim_dir")
    parser.add_argument("--step", type=int, required=True)
    parser.add_argument("--tol", type=float, default=1.0e-6)
    args = parser.parse_args(argv)

    path = os.path.join(args.sim_dir, "results", "elts", "stress", f"stress.step{args.step}")
    if not os.path.isfile(path):
        print(f"FAIL  missing {path} (did the run stop early?)")
        return 1
    # Stress rows: 11 22 33 23 13 12 (tensor components).
    stress = np.loadtxt(path, ndmin=2)
    other = np.max(np.abs(stress[:, [0, 1, 3, 4, 5]]), axis=1)
    ratio = other / np.abs(stress[:, 2])
    worst = float(np.max(ratio))
    n_fail = int(np.count_nonzero(ratio > args.tol))
    print(f"check_uniaxial_stress: step {args.step}, {len(stress)} elements, mean sigma_zz "
          f"{stress[:, 2].mean():.6g}, worst |other|/|sigma_zz| {worst:.3g} (tol {args.tol:g}), "
          f"{n_fail} failed")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
