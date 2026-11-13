#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check the lattice rotation rate output against the orientations (CHANGES.md C41).

Source : FEPX verification tests (ROADMAP step 6.5). Development tool (CLAUDE.md §6).

The lattice rotation rate of an increment (v2 write_res `rotrate`: the rotation vector of the
increment's lattice rotation divided by the time increment, in the sample frame) follows from
the orientations at its start and end: with rot_c2s the crystal-to-sample rotations,

    rotrate = rotvec(R_n R_(n-1)^T) / (t_n - t_(n-1)),

where rotvec is the rotation vector (axis times angle, radians) of a rotation matrix. When every
step has a single increment, the values written at step n must equal this finite difference of
the orientations written at steps n - 1 and n, for every element. The orientations must be
written as Rodrigues vectors in the passive convention (the regression meshes); their matrix g
is sample-to-crystal and rot_c2s = g^T (conventions.md §2).

Usage: check_rotrate.py SIM_DIR --times T1 T2 ... [--var rotrate] [--rtol R] [--min-rate M]
  --times     end time of every step 1, 2, ... (one increment per step), step 0 at t = 0
  --rtol      allowed difference, relative to the largest finite-difference rate of the step
              (default 2e-2: .sim files of v2 have 7 digits; the C41 error is 25 % or more)
  --min-rate  the largest rate over all steps must exceed it (the check is not vacuous)
Exit status 0 if the check passes, 1 otherwise.

Functions:
    read_rows        rows of a result file
    rot_c2s_from_rodrigues_passive  crystal-to-sample rotations of passive Rodrigues vectors
    rotation_vector  rotation vectors of rotation matrices
    main             command-line interface
"""

import argparse
import os
import sys

import numpy as np


def read_rows(sim_dir, var, step):
    """Return the rows of results/elts/<var>/<var>.step<step> as a 2D array (one row per
    element, also for one value per element)."""
    path = os.path.join(sim_dir, "results", "elts", var, f"{var}.step{step}")
    return np.loadtxt(path, ndmin=2)


def rot_c2s_from_rodrigues_passive(rod):
    """Return the crystal-to-sample rotations (n, 3, 3) of passive Rodrigues vectors (n, 3).

    The passive matrix is g = ((1 - r.r) I + 2 r r^T - 2 [r]x) / (1 + r.r); rot_c2s = g^T.
    """
    rot = np.empty((len(rod), 3, 3))
    for i, r in enumerate(rod):
        rr = r @ r
        cross = np.array([[0.0, -r[2], r[1]], [r[2], 0.0, -r[0]], [-r[1], r[0], 0.0]])
        rot[i] = ((1.0 - rr) * np.eye(3) + 2.0 * np.outer(r, r) + 2.0 * cross) / (1.0 + rr)
    return rot


def rotation_vector(rot):
    """Return the rotation vectors (n, 3) (axis times angle, radians) of rotations (n, 3, 3).

    The axial vector of the skew part, (m32 - m23, m13 - m31, m21 - m12) / 2, is sin(angle)
    times the axis, and (tr m - 1) / 2 = cos(angle): the angle is atan2 of the two (accurate
    for small angles; the regression rotations per increment are far from 180 degrees).
    """
    out = np.zeros((len(rot), 3))
    for i, m in enumerate(rot):
        axial = 0.5 * np.array([m[2, 1] - m[1, 2], m[0, 2] - m[2, 0], m[1, 0] - m[0, 1]])
        sin_angle = np.linalg.norm(axial)
        if sin_angle > 0.0:
            angle = np.arctan2(sin_angle, 0.5 * (np.trace(m) - 1.0))
            out[i] = axial / sin_angle * angle
    return out


def main(argv=None):
    """Command-line interface; returns the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("sim_dir")
    parser.add_argument("--times", type=float, nargs="+", required=True)
    parser.add_argument("--var", default="rotrate")
    parser.add_argument("--rtol", type=float, default=2.0e-2)
    parser.add_argument("--min-rate", type=float, default=0.0)
    args = parser.parse_args(argv)

    times = [0.0] + args.times
    rot_prev = rot_c2s_from_rodrigues_passive(read_rows(args.sim_dir, "ori", 0))
    ok = True
    largest = 0.0
    for step in range(1, len(times)):
        rot = rot_c2s_from_rodrigues_passive(read_rows(args.sim_dir, "ori", step))
        dt = times[step] - times[step - 1]
        # Finite difference of the orientations: the increment's rotation in the sample frame.
        fd = rotation_vector(np.einsum("eij,ekj->eik", rot, rot_prev)) / dt
        written = read_rows(args.sim_dir, args.var, step)
        scale = np.max(np.abs(fd))
        largest = max(largest, scale)
        diff = np.max(np.abs(written - fd))
        good = diff <= args.rtol * scale
        ok = ok and good
        print(f"check_rotrate: step {step}, dt {dt:.6g} s, largest finite-difference rate "
              f"{scale:.4e} 1/s, largest difference of {args.var} {diff:.3e} "
              f"({diff / scale if scale > 0 else 0.0:.3e} relative; rtol {args.rtol:.1e}) -> "
              f"{'ok' if good else 'FAIL'}")
        rot_prev = rot
    if not largest > args.min_rate:
        print(f"check_rotrate: largest rate {largest:.3e} not above --min-rate "
              f"{args.min_rate:.3e} -> FAIL")
        ok = False
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
