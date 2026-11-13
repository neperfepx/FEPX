#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check the initial face velocities of a triaxial step change (FEPX step 8.5, ROADMAP D86): the
first boundary-velocity iteration of step 2 starts from the face velocities of step 1's last
iteration times the ratio of the control component's strain rates.

Source : v2 driver_triaxcsr (`results%vel` and `curr_vel` times `vel_factor` at a step's start);
         dev/theory/addendum_triaxial_load_control.md §3.4 (revised, D86). Development tool
         (CLAUDE.md §6: Python 3 + NumPy for comparing outputs in tests).

The run's log (run.log next to SIM_DIR; step 9.7) prints at every boundary-velocity iteration
  "Load iteration <k>: face velocities <v_1> <v_2> <v_3> mm/s; ..."
and "Running step <s>..." before each step. Let v^(1) be the face velocities
of the last iteration printed in step 1, and v^(2) those of the first iteration of step 2; with
r = --ratio (the control's strain rate in step 2 over step 1):
  every load-controlled component j (--load-comps): v^(2)_j = r v^(1)_j;
  the control component c (--control): v^(2)_c = r v^(1)_c (it is at its strain rate times the
  length in both steps, the length unchanged between them).
Tolerance: --rtol of max_j |r v^(1)_j| (the log prints the face velocities with six significant
digits).

Usage:
  check_step_change_velocities.py SIM_DIR --control 11 --load-comps 22 33 --ratio 2 [--rtol R]
Exit status 0 when every check passes, 1 otherwise (message on standard output).
"""

import argparse
import os
import re
import sys

import numpy as np

from fepx_log import joined_lines


def face_velocities(log_text):
    """Face velocities of the first and the last boundary-velocity iteration of each step.

    Returns a dict step -> (first, last), each a NumPy array of the three face velocities [mm/s].
    """
    step = None
    found = {}
    for line in joined_lines(log_text):
        running = re.search(r"Running step (\d+)\.\.\.", line)
        if running:
            step = int(running.group(1))
            continue
        iteration = re.search(r"Load iteration \d+: face velocities +(\S+) +(\S+) +(\S+)", line)
        if iteration and step is not None:
            values = np.array([float(v) for v in iteration.groups()])
            first, _ = found.get(step, (values, values))
            found[step] = (first, values)
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("sim_dir")
    parser.add_argument("--control", type=int, required=True)
    parser.add_argument("--load-comps", type=int, nargs="+", required=True)
    parser.add_argument("--ratio", type=float, required=True)
    parser.add_argument("--rtol", type=float, default=2e-5)
    args = parser.parse_args()

    # The log of the run.
    log_path = os.path.join(os.path.dirname(os.path.abspath(args.sim_dir)), "run.log")
    with open(log_path, encoding="utf-8") as stream:
        found = face_velocities(stream.read())
    if 1 not in found or 2 not in found:
        print(f"FAIL  no boundary-velocity iterations of steps 1 and 2 in {log_path}")
        return 1

    # Step 2's first iteration against step 1's last one times the ratio.
    expected = args.ratio*found[1][1]
    actual = found[2][0]
    tol = args.rtol*np.max(np.abs(expected))
    failed = 0
    for comp in [args.control] + args.load_comps:
        j = comp//10 - 1
        diff = abs(actual[j] - expected[j])
        status = "ok  " if diff <= tol else "FAIL"
        failed += diff > tol
        print(f"{status}  component {comp}: step 2 starts at {actual[j]: .5e} mm/s, "
              f"expected {args.ratio:g} x {found[1][1][j]: .5e} = {expected[j]: .5e} "
              f"(|diff| {diff:.2e}, tolerance {tol:.2e})")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
