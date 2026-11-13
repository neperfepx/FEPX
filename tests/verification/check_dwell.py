#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check the dwells of a run at a constant load rate (FEPX step 8.2c): the load held, creep, and
v2's dwell increments.

Source : FEPX verification tests (ROADMAP step 8.2c, dev/theory/addendum_triaxial_load_control.md
         §3.3 and §5 V4; v2 driver_triaxclr's dwell time stepping). Development tool (CLAUDE.md §6:
         Python 3 + NumPy for comparing outputs in tests).

Theory (addendum §3.3): a dwell (load_rate 0 on the control component c) holds the value of c at
its start; its increments are 0.01 dtime at the start of a dwell episode (a dwell after a step that
is not a dwell), then min(max_strain_incr L0 / |V|, dtime), V the control face velocity of the
previous increment and L0 the original length along c, the last increment the remaining time (when
less than 0.001 dtime would remain after a full one, all of it).

Inputs: the forces file of the control face (results/forces/<c>1: step, increment, force, time;
an increment's dt is the difference of its time and the previous one) and the run log
(WORK_DIR/run.log, beside the .sim directory; step 9.7): the lines "Increment <n>: t = ..." and the
last "Load iteration ...: face velocities vx vy vz mm/s; ..." of each increment (its accepted face
velocities).

Checks, for each dwell step in --dwells (with its dtime and end time):
  1. every increment's control load within --load-tol of the load at the dwell's start;
  2. the control face moves on in the loading direction in every increment (creep: V s > 0, s the
     sign of the held load);
  3. every increment's dt is v2's rule (relative 2e-5: the log writes V with 6 significant
     digits), the dwell ending at its time (1e-9 s, from the forces file).

Usage:
  check_dwell.py SIM_DIR --component 33 --dwells 2:0.5:4.2 3:0.5:5.2 --episode-starts 2
      --max-strain-incr 2e-4 --length0 1.0 --load-tol 0.1
  (each dwell as step:dtime:end time; --episode-starts the dwell steps that start an episode)
Exit status 0 when every check passes, 1 otherwise (message on standard output).
"""

import argparse
import os
import re
import sys

import numpy as np

from fepx_log import joined_lines


def read_log(path, comp):
    """Increment -> (step, the control face velocity of its last iteration)."""
    increments = {}
    current = None
    step = None
    with open(path, encoding="utf-8", errors="replace") as log:
        for line in joined_lines(log.read()):
            m = re.search(r"Running step (\d+)\.\.\.", line)
            if m:
                step = int(m.group(1))
                continue
            m = re.search(r"Increment (\d+): t = ", line)
            if m:
                current = int(m.group(1))
                increments[current] = [step, None]
                continue
            m = re.search(r"Load iteration \d+: face velocities +(\S+) +(\S+) +(\S+) mm/s", line)
            if m and current is not None:
                increments[current][1] = float(m.group(comp))
    return increments


def main(argv=None):
    """Command-line interface; returns the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("sim_dir")
    parser.add_argument("--component", type=int, choices=(11, 22, 33), required=True)
    parser.add_argument("--dwells", nargs="+", required=True, help="step:dtime:end_time")
    parser.add_argument("--episode-starts", type=int, nargs="*", default=[])
    parser.add_argument("--max-strain-incr", type=float, required=True)
    parser.add_argument("--length0", type=float, required=True)
    parser.add_argument("--load-tol", type=float, required=True)
    args = parser.parse_args(argv)

    comp = args.component // 11
    face = "xyz"[comp - 1] + "1"
    table = np.loadtxt(os.path.join(args.sim_dir, "results", "forces", face), comments="%",
                       ndmin=2)
    steps, incrs, force, time = table[:, 0].astype(int), table[:, 1].astype(int), \
        table[:, 1 + comp], table[:, 6]
    log = read_log(os.path.join(args.sim_dir, os.pardir, "run.log"), comp)

    failures = []
    for spec in args.dwells:
        k, dtime, t_end = spec.split(":")
        k, dtime, t_end = int(k), float(dtime), float(t_end)
        idx = np.nonzero(steps == k)[0]
        start = idx[0] - 1
        f0 = force[start]
        # 1. The load held.
        drift = np.max(np.abs(force[idx] - f0))
        if drift > args.load_tol * (1 + 1e-9):
            failures.append(f"step {k}: load off its start value {f0:.6e} by {drift:.3e} N")
        # 2. Creep, 3. the increments.
        for j, i in enumerate(idx):
            n = incrs[i]
            _, vel = log[n]
            dt = time[i] - time[i - 1]
            if not vel * np.sign(f0) > 0.0:
                failures.append(f"step {k}, increment {n}: face velocity {vel:.3e}: no creep")
            remaining = t_end - time[i - 1]
            if j == 0 and k in args.episode_starts:
                expect = 0.01 * dtime
            else:
                v_prev = log[incrs[i - 1]][1]
                expect = min(args.max_strain_incr * args.length0 / abs(v_prev), dtime)
            if not remaining - expect > 0.001 * dtime:
                expect = remaining
            if abs(dt - expect) > 2e-5 * expect:
                failures.append(f"step {k}, increment {n}: dt {dt:.9e} s, v2's rule {expect:.9e}")
        if abs(time[idx[-1]] - t_end) > 1e-9:
            failures.append(f"step {k}: ends at {time[idx[-1]]:.12e} s, not {t_end}")
        print(f"step {k}: {idx.size} dwell increments, load drift {drift:.3e} N, "
              f"dt {[round(time[i] - time[i - 1], 6) for i in idx]}")

    for message in failures:
        print("FAIL: " + message)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
