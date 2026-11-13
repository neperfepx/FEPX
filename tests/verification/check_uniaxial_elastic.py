#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check a single-crystal uniaxial-stress elastic run against the analytic anisotropic compliance.

Source : FEPX verification tests (ROADMAP step 0.4, CHANGES.md C1). Development tool
         (CLAUDE.md §6: Python 3 + NumPy for comparing outputs in tests).

Theory (linear anisotropic elasticity, small strain): with uniaxial stress sigma_kk along sample axis
k (--load-dir 1, 2 or 3 for x, y, z; default 3), the crystal-frame stress is sigma_c = R^T sigma_s R, the crystal-frame strain is
eps_c = [C_c]^-1 sigma_c, and the sample-frame strain is eps_s = R eps_c R^T, where R is the
crystal-to-sample rotation (v_s = R v_c). [C_c] is the hexagonal/tetragonal stiffness with the
c-axis along crystal axis 3 and Voigt constants c11 c12 c13 c33 c44 c66.

FEPX v2's deviatoric/volumetric split constrains c33 = c11 + c12 - c13 (and, for hcp,
c66 = (c11 - c12)/2); these are the defaults when --c33/--c66 are not given.

Checks, for every element of step --step (which must exist: v2 can stop early with exit status 0,
CHANGES.md C10):
  1. the stress is uniaxial: every component other than sigma_kk is below --stress-tol * |sigma_kk|;
  2. every strain component agrees with the analytic strain for that element's sigma_kk:
     |eps - eps_analytic| <= --rtol * max|eps_analytic|.

Usage:
  check_uniaxial_elastic.py SIM_DIR --step N --c11 .. --c12 .. --c13 .. --c44 .. [--c33 ..]
                            [--c66 ..] --axis X Y Z --angle DEG [--load-dir K] [--rtol R]
                            [--stress-tol T] [--strain-var strain|strain_el]
Exit status 0 if all checks pass, 1 otherwise.

Functions:
    stiffness_hex      3x3x3x3 stiffness of a crystal with hexagonal/tetragonal symmetry
    rotation_axis_angle  rotation matrix from an axis and an angle (right-hand rule)
    analytic_strain    sample-frame strain (6 tensor components, order 11 22 33 23 13 12)
    read_step          rows of an element result variable at one step
    main               command-line interface
"""

import argparse
import os
import sys

import numpy as np

# Tensor component order of FEPX output: 11 22 33 23 13 12.
ORDER = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
# Mandel weights for that order.
WEIGHTS = np.array([1.0, 1.0, 1.0, np.sqrt(2.0), np.sqrt(2.0), np.sqrt(2.0)])


def stiffness_hex(c11, c12, c13, c33, c44, c66):
    """Return the 3x3x3x3 stiffness (crystal frame, c-axis = axis 3) from Voigt constants."""
    stiff = np.zeros((3, 3, 3, 3))

    def put(i, j, k, l, value):
        # Fill all minor and major symmetric partners.
        for a, b in ((i, j), (j, i)):
            for c, d in ((k, l), (l, k)):
                stiff[a, b, c, d] = value
                stiff[c, d, a, b] = value

    put(0, 0, 0, 0, c11)
    put(1, 1, 1, 1, c11)
    put(2, 2, 2, 2, c33)
    put(0, 0, 1, 1, c12)
    put(0, 0, 2, 2, c13)
    put(1, 1, 2, 2, c13)
    put(1, 2, 1, 2, c44)
    put(0, 2, 0, 2, c44)
    put(0, 1, 0, 1, c66)
    return stiff


def rotation_axis_angle(axis, angle_deg):
    """Rotation matrix for a right-handed rotation of angle_deg about axis (Rodrigues formula)."""
    unit = np.asarray(axis, dtype=float)
    unit = unit / np.linalg.norm(unit)
    theta = np.radians(angle_deg)
    skew = np.array([[0.0, -unit[2], unit[1]], [unit[2], 0.0, -unit[0]], [-unit[1], unit[0], 0.0]])
    return np.eye(3) + np.sin(theta) * skew + (1.0 - np.cos(theta)) * (skew @ skew)


def analytic_strain(stiff_c, rot_c2s, sigma_kk, load_dir=3):
    """Sample-frame strain (6 tensor components) under uniaxial stress sigma_kk along axis k."""
    # Mandel 6x6 stiffness in the crystal frame.
    mandel = np.array([[stiff_c[i, j, k, l] * WEIGHTS[a] * WEIGHTS[b]
                        for b, (k, l) in enumerate(ORDER)] for a, (i, j) in enumerate(ORDER)])
    # Stress: sample -> crystal frame.
    sigma_s = np.zeros((3, 3))
    sigma_s[load_dir - 1, load_dir - 1] = sigma_kk
    sigma_c = rot_c2s.T @ sigma_s @ rot_c2s
    # Crystal-frame strain from the compliance.
    sig_vec = np.array([sigma_c[i, j] * WEIGHTS[a] for a, (i, j) in enumerate(ORDER)])
    eps_vec = np.linalg.solve(mandel, sig_vec) / WEIGHTS
    eps_c = np.zeros((3, 3))
    for a, (i, j) in enumerate(ORDER):
        eps_c[i, j] = eps_c[j, i] = eps_vec[a]
    # Strain: crystal -> sample frame.
    eps_s = rot_c2s @ eps_c @ rot_c2s.T
    return np.array([eps_s[i, j] for i, j in ORDER])


def read_step(sim_dir, var, step):
    """Return the rows of results/elts/<var>/<var>.step<step> (error if missing)."""
    path = os.path.join(sim_dir, "results", "elts", var, f"{var}.step{step}")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"missing {path} (did the run stop early?)")
    return np.loadtxt(path, ndmin=2)


def main(argv=None):
    """Command-line interface; returns the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("sim_dir")
    parser.add_argument("--step", type=int, required=True, help="step to check (must exist)")
    for name in ("c11", "c12", "c13", "c44"):
        parser.add_argument(f"--{name}", type=float, required=True)
    parser.add_argument("--c33", type=float, help="default: c11 + c12 - c13 (v2 constraint)")
    parser.add_argument("--c66", type=float, help="default: (c11 - c12)/2 (hexagonal)")
    parser.add_argument("--axis", type=float, nargs=3, required=True)
    parser.add_argument("--angle", type=float, required=True, help="degrees")
    parser.add_argument("--load-dir", type=int, choices=(1, 2, 3), default=3,
                        help="loading axis: 1 x, 2 y, 3 z (default 3)")
    parser.add_argument("--rtol", type=float, default=1.0e-3)
    parser.add_argument("--stress-tol", type=float, default=1.0e-3)
    parser.add_argument("--strain-var", default="strain",
                        help="strain variable to check: strain (v2) or strain_el (FEPX, step 4.4:"
                             " one elastic increment from zero, where it equals the total strain)")
    args = parser.parse_args(argv)

    # Elastic constants (v2 constraints by default) and rotation.
    c33 = args.c33 if args.c33 is not None else args.c11 + args.c12 - args.c13
    c66 = args.c66 if args.c66 is not None else (args.c11 - args.c12) / 2.0
    stiff_c = stiffness_hex(args.c11, args.c12, args.c13, c33, args.c44, c66)
    rot_c2s = rotation_axis_angle(args.axis, args.angle)

    # Results of the requested step.
    step = args.step
    try:
        stress = read_step(args.sim_dir, "stress", step)
        strain = read_step(args.sim_dir, args.strain_var, step)
    except FileNotFoundError as err:
        print(f"FAIL  {err}")
        return 1

    n_fail = 0
    worst = 0.0
    for i_elt, (sig, eps) in enumerate(zip(stress, strain), start=1):
        k = args.load_dir - 1
        sigma_kk = sig[k]
        # 1. Uniaxial stress.
        other = np.max(np.abs(np.delete(sig, k)))
        if other > args.stress_tol * abs(sigma_kk):
            n_fail += 1
            print(f"FAIL  elt {i_elt}: stress not uniaxial: {sig}")
            continue
        # 2. Strain against the analytic compliance.
        eps_ref = analytic_strain(stiff_c, rot_c2s, sigma_kk, args.load_dir)
        err = np.max(np.abs(eps - eps_ref)) / np.max(np.abs(eps_ref))
        worst = max(worst, err)
        if err > args.rtol:
            n_fail += 1
            print(f"FAIL  elt {i_elt}: strain {eps} vs analytic {eps_ref} (rel. err {err:.3g})")
    print(f"check_uniaxial_elastic: step {step}, {len(stress)} elements, {n_fail} failed, "
          f"worst strain error {worst:.3g} (rtol {args.rtol:g}), c33 {c33:g}, c66 {c66:g}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
