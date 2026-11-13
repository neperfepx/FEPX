#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check the steady power-law flow stress of a single crystal under uniaxial stress.

Source : FEPX verification tests (ROADMAP step 6.2). Development tool (CLAUDE.md §6).

A single crystal loaded along z (uniaxial stress) whose plastic flow is carried by N systems with
the same Schmid factor S and the same strength crss reaches, once the elastic strain rate vanishes,
d_zz = N S gammadot on each system, so that the Kirchhoff stress is

    tau_zz = crss (d_zz / (N S gammadot_0))^m / S,

with d_zz the true (logarithmic) strain rate, rate / (1 + strain) for an engineering strain rate
at an engineering strain. The Cauchy stress written by FEPX is tau / beta, beta = det(I + e^e),
computed in each element from its elastic strain (Mandel, 11 22 33 23 13 12). Every element's
sigma_zz at step --step must equal that value within --tol (relative).

Usage: check_flow_stress.py SIM_DIR --step N --crss C --m M --gammadot0 G0 --rate R --strain E
                            --schmid S --num-systems N [--tol T]
Exit status 0 if the check passes, 1 otherwise.
"""

import argparse
import math
import os
import sys

import numpy as np


def main(argv=None):
    """Command-line interface; returns the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("sim_dir")
    for name in ("step", "num-systems"):
        parser.add_argument(f"--{name}", type=int, required=True)
    for name in ("crss", "m", "gammadot0", "rate", "strain", "schmid"):
        parser.add_argument(f"--{name}", type=float, required=True)
    parser.add_argument("--tol", type=float, default=1.0e-5)
    args = parser.parse_args(argv)

    base = os.path.join(args.sim_dir, "results", "elts")
    stress = np.atleast_2d(np.loadtxt(os.path.join(base, "stress", f"stress.step{args.step}")))
    strain_el = np.atleast_2d(np.loadtxt(os.path.join(base, "strain_el",
                                                      f"strain_el.step{args.step}")))

    # Steady Kirchhoff stress at the true strain rate.
    d_zz = args.rate / (1.0 + args.strain)
    tau = args.crss * (d_zz / (args.num_systems * args.schmid * args.gammadot0)) ** args.m \
        / args.schmid

    s2 = math.sqrt(2.0)
    worst = 0.0
    for sig, e in zip(stress, strain_el):
        mat = np.array([[e[0], e[5] / s2, e[4] / s2],
                        [e[5] / s2, e[1], e[3] / s2],
                        [e[4] / s2, e[3] / s2, e[2]]])
        expected = tau / np.linalg.det(np.eye(3) + mat)
        worst = max(worst, abs(sig[2] / expected - 1.0))
    ok = worst <= args.tol
    print(f"{'ok' if ok else 'FAIL'}: {len(stress)} elements, largest relative difference of "
          f"sigma_zz from the steady flow stress {worst:.3e} (tolerance {args.tol:.1e})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
