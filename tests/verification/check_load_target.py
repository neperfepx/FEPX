#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check the load and stress targets of an elastic run against v2's step-completion rule.

Source : FEPX verification tests (ROADMAP step 8.1, dev/theory/addendum_load_targets.md §2,
         §5 V1-V2). Development tool (CLAUDE.md §6: Python 3 + NumPy for comparing outputs in
         tests).

Theory (addendum §2, v2 loading_timeincr / loading_load_isinrange): a load (force on the face
of the control component, column F<c> of results/forces/<face>) or stress (that force over the
face's current area) target step ends after the first increment whose value v satisfies
(v - v*) s >= -tol, s = sign(v* - v_start); every increment is at most dt_max and at least
dt_min; when the value extrapolated over dt_max with the last increment's rate reaches the
target, the next increment is f (v* - v)/vdot. For a linear response (elasticity, small
strains) the extrapolation is exact, so the step's last increment, when aimed (shorter than
dt_max), ends at v_end = v_prev + f (v* - v_prev): the overshoot is (f - 1)(v* - v_prev).

Since step 8.2a the tolerance tol is load_tol [N] or stress_tol [MPa] (--load-tol,
--stress-tol), and an increment that passes the target by more than tol is redone (landing rule,
addendum_triaxial_load_control.md §3.6; checked by check_landing.py): this check needs a case
whose tolerances exceed the overshoots, so that no increment is redone.

Checks, for every step of --targets (kind:value, in step order):
  1. the step ends at its target (within the tolerance short of it, or past it by at most the
     tolerance) and no earlier increment of the step reached it;
  2. every increment of the step lies in [dt_min, dt_max] (to 1e-9 relative);
  3. an aimed last increment overshoots by (f - 1)(v* - v_prev), to --rtol |v* - v_prev| (the
     residual nonlinearity of the small-strain geometry update).

Exit status 0 when every check passes, 1 otherwise (message on standard output).
"""

import argparse
import os
import sys

import numpy as np


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("sim_dir")
    parser.add_argument("--targets", nargs="+", required=True,
                        help="kind:value per step, kind load [N] or stress [MPa]")
    parser.add_argument("--face", default="z1", help="control face set (default z1)")
    parser.add_argument("--comp", type=int, choices=(1, 2, 3), default=3,
                        help="control component (default 3)")
    parser.add_argument("--factor", type=float, required=True, help="dtime_factor")
    parser.add_argument("--dt-max", type=float, nargs="+", required=True,
                        help="dtime of each step (one value: every step) [s]")
    parser.add_argument("--dt-min", type=float, required=True, help="dtime_min of every step [s]")
    parser.add_argument("--load-tol", type=float, default=0.1, help="load_tol [N]")
    parser.add_argument("--stress-tol", type=float, default=0.1, help="stress_tol [MPa]")
    parser.add_argument("--rtol", type=float, default=1.0e-3)
    args = parser.parse_args(argv)

    # The force table: step incr Fx Fy Fz area time, one row per increment (row 0: step 0).
    table = np.loadtxt(os.path.join(args.sim_dir, "results", "forces", args.face),
                       comments="%", ndmin=2)
    steps = table[:, 0].astype(int)
    force = table[:, 1 + args.comp]
    area = table[:, 5]
    time = table[:, 6]

    failures = []
    for k, spec in enumerate(args.targets, start=1):
        kind, value = spec.split(":")
        target = float(value)
        rows = np.nonzero(steps == k)[0]
        if rows.size == 0:
            failures.append(f"step {k}: no increments in the force table")
            continue
        values = force / area if kind == "stress" else force
        v_start = values[rows[0] - 1]
        sign = np.sign(target - v_start)
        tol = args.stress_tol if kind == "stress" else args.load_tol

        # 1. Completion at the last increment (not passed by more than tol), not before.
        reached = (values[rows] - target) * sign >= -tol
        if abs(values[rows[-1]] - target) > tol:
            failures.append(f"step {k}: ends at {values[rows[-1]]:.6e}, more than {tol} from "
                            f"the target {target:.6e}")
        if not reached[-1]:
            failures.append(f"step {k}: ends at {values[rows[-1]]:.6e}, target {target:.6e} "
                            "not reached")
        if np.any(reached[:-1]):
            failures.append(f"step {k}: target reached before the last increment")

        # 2. Increment bounds.
        dt_max = args.dt_max[min(k, len(args.dt_max)) - 1]
        dts = time[rows] - time[rows - 1]
        if np.any(dts > dt_max * (1 + 1e-9)) or np.any(dts < args.dt_min * (1 - 1e-9)):
            failures.append(f"step {k}: increments {dts} outside [{args.dt_min}, {dt_max}]")

        # 3. Overshoot of an aimed last increment (needs a previous increment in the step for
        #    the rate).
        v_prev = values[rows[-1] - 1]
        if rows.size >= 2 and dts[-1] < dt_max * (1 - 1e-9):
            expected = v_prev + args.factor * (target - v_prev)
            error = abs(values[rows[-1]] - expected)
            scale = abs(target - v_prev)
            status = "ok" if error <= args.rtol * scale else "FAIL"
            print(f"step {k} ({kind} {target:g}): {rows.size} increments, end "
                  f"{values[rows[-1]]:.8e}, expected {expected:.8e}, error/|v*-v_prev| "
                  f"{error / scale:.2e} [{status}]")
            if status == "FAIL":
                failures.append(f"step {k}: overshoot not (f - 1)(v* - v_prev)")
        else:
            failures.append(f"step {k}: the last increment is not an aimed one (case set-up)")

    for message in failures:
        print("FAIL: " + message)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
