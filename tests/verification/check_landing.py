#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check the landing of load and stress targets (FEPX step 8.2a) in a run's forces and log.

Source : FEPX verification tests (ROADMAP step 8.2a, dev/theory/addendum_triaxial_load_control.md
         §3.6 and §5 V7). Development tool (CLAUDE.md §6: Python 3 + NumPy for comparing outputs
         in tests).

Theory (addendum §3.6): a load (force on the face of the control component, column F<c> of
results/forces/<face>) or stress (that force over the face's current area) target step ends at
the first increment whose value v is within the tolerance tol (load_tol [N], stress_tol [MPa])
of the target v*; an increment that passes the target by more than tol is discarded and redone
with the length given by false position on the bracket (0, v_start of the increment),
(dt, v_passed). Only the accepted attempt is written. When the value is linear in the increment
length (elasticity, small strains), the first attempt lands on the target up to the residual
nonlinearity of the geometry update and the solver tolerance.

Checks, for every step of --targets (kind:value, in step order):
  1. the step ends within tol of its target, and every earlier increment of the step is short
     of it by more than tol;
  2. the log (run.log next to the .sim directory) reports, for the step's last increment, the
     number of landing attempts given by --attempts (0: the increment was not redone, and no
     landing is reported);
  3. a landed increment ends at the target to --land-rtol |v_end - v_prev| (v_prev the value
     before the increment);
  4. every increment of the step is at most dt_max (to 1e-9 relative).

Exit status 0 when every check passes, 1 otherwise (message on standard output).
"""

import argparse
import os
import re
import sys

import numpy as np

from fepx_log import joined_lines


def landed_attempts(log_path):
    """Landing attempts per increment from the log: {increment: attempts}.

    The log (step 9.7) has "- Increment <n>: t = ..." and, under it, "> Landed at <value> <unit>
    in <k> attempt(s) (dt = ...)".
    """
    increment = re.compile(r"Increment (\d+): t = ")
    landed = re.compile(r"Landed at .* in (\d+) attempts? ")
    attempts = {}
    current = None
    with open(log_path, errors="replace") as log:
        for line in joined_lines(log.read()):
            match = increment.search(line)
            if match:
                current = int(match.group(1))
                continue
            match = landed.search(line)
            if match and current is not None:
                attempts[current] = int(match.group(1))
    return attempts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("sim_dir")
    parser.add_argument("--targets", nargs="+", required=True,
                        help="kind:value per step, kind load [N] or stress [MPa]")
    parser.add_argument("--attempts", type=int, nargs="+", required=True,
                        help="landing attempts of each step's last increment (0: none)")
    parser.add_argument("--face", default="z1", help="control face set (default z1)")
    parser.add_argument("--comp", type=int, choices=(1, 2, 3), default=3,
                        help="control component (default 3)")
    parser.add_argument("--dt-max", type=float, nargs="+", required=True,
                        help="dtime of each step (one value: every step) [s]")
    parser.add_argument("--load-tol", type=float, default=0.1, help="load_tol [N]")
    parser.add_argument("--stress-tol", type=float, default=0.1, help="stress_tol [MPa]")
    parser.add_argument("--land-rtol", type=float, default=1.0e-4)
    args = parser.parse_args(argv)
    if len(args.attempts) != len(args.targets):
        parser.error("give one --attempts value per target")

    # The force table: step incr Fx Fy Fz area time, one row per increment (row 0: step 0).
    table = np.loadtxt(os.path.join(args.sim_dir, "results", "forces", args.face),
                       comments="%", ndmin=2)
    steps = table[:, 0].astype(int)
    incrs = table[:, 1].astype(int)
    force = table[:, 1 + args.comp]
    area = table[:, 5]
    time = table[:, 6]
    attempts = landed_attempts(os.path.join(os.path.dirname(os.path.abspath(args.sim_dir)),
                                            "run.log"))

    failures = []
    for k, spec in enumerate(args.targets, start=1):
        kind, value = spec.split(":")
        target = float(value)
        tol = args.stress_tol if kind == "stress" else args.load_tol
        rows = np.nonzero(steps == k)[0]
        if rows.size == 0:
            failures.append(f"step {k}: no increments in the force table")
            continue
        values = force / area if kind == "stress" else force
        sign = np.sign(target - values[rows[0] - 1])

        # 1. Within tol at the end, short by more than tol before.
        if abs(values[rows[-1]] - target) > tol:
            failures.append(f"step {k}: ends at {values[rows[-1]]:.6e}, more than {tol} from "
                            f"the target {target:.6e}")
        if np.any((values[rows[:-1]] - target) * sign >= -tol):
            failures.append(f"step {k}: an increment before the last is within {tol} of the "
                            "target or past it")

        # 2. Landing attempts of the last increment, from the log.
        found = attempts.get(int(incrs[rows[-1]]), 0)
        if found != args.attempts[k - 1]:
            failures.append(f"step {k}: increment {incrs[rows[-1]]} landed in {found} "
                            f"attempt(s), expected {args.attempts[k - 1]}")

        # 3. A landed increment ends at the target (linear response).
        v_prev = values[rows[-1] - 1]
        error = abs(values[rows[-1]] - target)
        scale = abs(values[rows[-1]] - v_prev)
        status = "ok"
        if found > 0 and error > args.land_rtol * scale:
            status = "FAIL"
            failures.append(f"step {k}: landed {error:.3e} from the target, more than "
                            f"{args.land_rtol} x {scale:.3e}")
        print(f"step {k} ({kind} {target:g}): {rows.size} increments, end "
              f"{values[rows[-1]]:.8e}, {found} landing attempt(s), |v - v*|/|v - v_prev| "
              f"{error / scale:.2e} [{status}]")

        # 4. Increments at most dt_max.
        dt_max = args.dt_max[min(k, len(args.dt_max)) - 1]
        dts = time[rows] - time[rows - 1]
        if np.any(dts > dt_max * (1 + 1e-9)):
            failures.append(f"step {k}: increments {dts} above dt_max {dt_max}")

    for message in failures:
        print("FAIL: " + message)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
