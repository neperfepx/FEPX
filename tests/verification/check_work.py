#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check the work output of a homogeneous uniaxial run (step 6.5).

Source : FEPX verification tests (ROADMAP step 6.5). Development tool (CLAUDE.md §6).

Under uniaxial stress along z with the top face moving at rate * L0 (uniaxial_minimal, constant
engineering strain rate), the deformation is homogeneous and the deformation rate along z at
time t is d_zz = rate / (1 + rate t) (the true strain rate); the stress is uniaxial, so the work
rate sigma : d is sigma_zz d_zz. The work written at step n (one increment per step from
--first-step on; step 1 may contain the initialization increment) must follow the trapezoidal
rule of FEPX and v2 (finalize_res_totalwork):

    work_n = work_(n-1) + dt/2 (sigma_zz,n d_zz(t_n) + sigma_zz,(n-1) d_zz(t_(n-1))),

in every element, with sigma_zz from the stress output (the off-axis stress components are
checked to be below 1e-6 sigma_zz, so the neglected terms of sigma : d are negligible).

Usage: check_work.py SIM_DIR --times T1 T2 ... --rate R [--first-step K] [--tol T]
  --tol  allowed difference relative to the largest work increment of the step (default 1e-6:
         the deformation is homogeneous and the stress uniaxial only to the global solver
         tolerance; the neglected off-axis terms of sigma : d are of the size of the off-axis
         stress, about 1e-7 sigma_zz)
Exit status 0 if the check passes, 1 otherwise.
"""

import argparse
import sys

import numpy as np

from check_rotrate import read_rows


def main(argv=None):
    """Command-line interface; returns the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("sim_dir")
    parser.add_argument("--times", type=float, nargs="+", required=True)
    parser.add_argument("--rate", type=float, required=True)
    parser.add_argument("--first-step", type=int, default=1)
    parser.add_argument("--tol", type=float, default=1.0e-6)
    args = parser.parse_args(argv)

    times = [0.0] + args.times
    ok = True
    for step in range(args.first_step, len(times)):
        dt = times[step] - times[step - 1]
        stress = read_rows(args.sim_dir, "stress", step)
        stress_prev = read_rows(args.sim_dir, "stress", step - 1)

        # Uniaxial stress (the work rate is then sigma_zz d_zz).
        off_axis = np.max(np.abs(stress[:, [0, 1, 3, 4, 5]])) / np.max(np.abs(stress[:, 2]))
        good = off_axis < 1.0e-6
        print(f"check_work: step {step}, off-axis stress / sigma_zz {off_axis:.3e} -> "
              f"{'ok' if good else 'FAIL'}")
        ok = ok and good

        # Trapezoidal rule with the work rates at the start and end of the increment.
        rate_end = stress[:, 2] * args.rate / (1.0 + args.rate * times[step])
        rate_start = stress_prev[:, 2] * args.rate / (1.0 + args.rate * times[step - 1])
        work = read_rows(args.sim_dir, "work", step)[:, 0]
        work_prev = read_rows(args.sim_dir, "work", step - 1)[:, 0]
        expected = work_prev + 0.5 * dt * (rate_start + rate_end)
        diff = np.max(np.abs(work - expected))
        scale = np.max(np.abs(work - work_prev))
        good = diff <= args.tol * scale
        print(f"check_work: step {step}, dt {dt:.6g} s, largest work increment {scale:.4e} MPa, "
              f"largest difference {diff:.3e} ({diff / scale:.3e} relative; tol "
              f"{args.tol:.1e}) -> {'ok' if good else 'FAIL'}")
        ok = ok and good
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
