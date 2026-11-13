#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check a single-crystal triaxial strain-controlled elastic run against the homogeneous solution.

Source : FEPX verification tests (ROADMAP step 5.3, D45). Development tool (CLAUDE.md §6:
         Python 3 + NumPy for comparing outputs in tests).

Theory (linear anisotropic elasticity; recoupled_formulation.md R3-R4): the triaxial preset
prescribes the normal velocity of faces x0, y0, z0 (zero) and of the faces x1, y1, z1 of the
prescribed components, and leaves the shear tractions free. When the crystal's sample-frame
stiffness does not couple normal and shear components (a crystal axis along each sample axis,
e.g. hexagonal with the c-axis along a sample axis), the solution is homogeneous: zero shear
strain and stress, zero normal stress on the free components, and

    tau = [C_s] {e},  sigma = tau / beta,  beta = det(I + e),

with the free normal strains e_F = -[C_FF]^-1 [C_FP] e_P (P: prescribed normal components).
[C_s] = [Q6(R)][C_c][Q6(R)]^T, R the crystal-to-sample rotation. (With normal-shear coupling
the faces, held plane with normal motion only, would need shear tractions: no homogeneous
solution, so such crystals are not checked here.)

Checks, for every element of step --step (which must exist), with the computed elastic strain e
(--strain-var, default strain_el) and stress sigma:
  1. shear strains below --rtol * max |e|, shear stresses and the normal stresses of the free
     components below --stress-tol * max |sigma|;
  2. sigma = [C_s]{e}/beta to --rtol * max |sigma|;
  3. the free normal strains equal -[C_FF]^-1 [C_FP] e_P to --rtol * max |e|;
  4. with --expect, the prescribed normal strains equal the expected engineering strains to
     --target-rtol * max |expected| (the elastic strain accumulates d dt, which differs from
     the engineering strain by about e/2 relative).

Usage:
  check_triaxial_elastic.py SIM_DIR --step N --c11 .. --c12 .. --c13 .. --c33 .. --c44 ..
                            [--c66 ..] --axis X Y Z --angle DEG --prescribed 11 [22] [33]
                            [--expect=E11,E22,E33] [--rtol R] [--stress-tol T]
                            [--target-rtol T]
Exit status 0 if all checks pass, 1 otherwise.

Functions:
    mandel_stiffness  sample-frame Mandel stiffness of the crystal
    main              command-line interface
"""

import argparse
import sys

import numpy as np

from check_uniaxial_elastic import ORDER, WEIGHTS, read_step, rotation_axis_angle, stiffness_hex


def mandel_stiffness(stiff_c, rot_c2s):
    """Return the sample-frame 6x6 Mandel stiffness (order 11 22 33 23 13 12)."""
    # Rotate the 4th-order tensor, then build the Mandel matrix.
    stiff_s = np.einsum("ia,jb,kc,ld,abcd->ijkl", rot_c2s, rot_c2s, rot_c2s, rot_c2s, stiff_c)
    return np.array([[stiff_s[i, j, k, l] * WEIGHTS[a] * WEIGHTS[b]
                      for b, (k, l) in enumerate(ORDER)] for a, (i, j) in enumerate(ORDER)])


def main(argv=None):
    """Command-line interface; returns the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("sim_dir")
    parser.add_argument("--step", type=int, required=True)
    for name in ("c11", "c12", "c13", "c33", "c44"):
        parser.add_argument(f"--{name}", type=float, required=True)
    parser.add_argument("--c66", type=float, help="default: (c11 - c12)/2 (hexagonal)")
    parser.add_argument("--axis", type=float, nargs=3, required=True)
    parser.add_argument("--angle", type=float, required=True, help="degrees")
    parser.add_argument("--prescribed", type=int, nargs="+", required=True,
                        choices=(11, 22, 33), help="prescribed normal components")
    parser.add_argument("--expect", type=lambda text: [float(v) for v in text.split(",")],
                        help="expected engineering strains 'e11,e22,e33' (prescribed ones used;"
                             " one comma-separated value written --expect=..., since argparse"
                             " before Python 3.14 takes a leading '-2e-4' for an option)")
    parser.add_argument("--rtol", type=float, default=1.0e-6)
    parser.add_argument("--stress-tol", type=float, default=1.0e-6)
    parser.add_argument("--target-rtol", type=float, default=3.0e-3)
    parser.add_argument("--strain-var", default="strain_el")
    args = parser.parse_args(argv)

    # Sample-frame stiffness; prescribed and free normal components (indices 0..2).
    c66 = args.c66 if args.c66 is not None else (args.c11 - args.c12) / 2.0
    stiff_c = stiffness_hex(args.c11, args.c12, args.c13, args.c33, args.c44, c66)
    rot_c2s = rotation_axis_angle(args.axis, args.angle)
    stiff_s = mandel_stiffness(stiff_c, rot_c2s)
    pres = sorted(c // 11 - 1 for c in args.prescribed)
    free = [k for k in range(3) if k not in pres]
    if np.max(np.abs(stiff_s[:3, 3:])) > 1.0e-9 * np.max(np.abs(stiff_s)):
        print("FAIL  the crystal couples normal and shear components in the sample frame: "
              "no homogeneous solution to check")
        return 1

    # Results of the requested step (tensor components, order 11 22 33 23 13 12).
    try:
        stress = read_step(args.sim_dir, "stress", args.step)
        strain = read_step(args.sim_dir, args.strain_var, args.step)
    except FileNotFoundError as err:
        print(f"FAIL  {err}")
        return 1

    n_fail = 0
    worst = {"shear strain": 0.0, "shear/free stress": 0.0, "hooke": 0.0, "free strain": 0.0,
             "target": 0.0}
    for i_elt, (sig, eps) in enumerate(zip(stress, strain), start=1):
        eps_scale = np.max(np.abs(eps[:3]))
        sig_scale = np.max(np.abs(sig))
        errors = {}

        # 1. Shear strains, shear stresses and free normal stresses vanish.
        errors["shear strain"] = np.max(np.abs(eps[3:])) / eps_scale
        errors["shear/free stress"] = (np.max(np.abs(np.concatenate((sig[3:], sig[free]))))
                                       / sig_scale)
        # 2. Hooke's law with beta (Mandel vectors from the tensor components).
        eps_mat = np.zeros((3, 3))
        for a, (i, j) in enumerate(ORDER):
            eps_mat[i, j] = eps_mat[j, i] = eps[a]
        beta = np.linalg.det(np.eye(3) + eps_mat)
        sig_ref = (stiff_s @ (eps * WEIGHTS)) / beta / WEIGHTS
        errors["hooke"] = np.max(np.abs(sig - sig_ref)) / sig_scale
        # 3. Free normal strains from the prescribed ones.
        if free:
            eps_free = -np.linalg.solve(stiff_s[np.ix_(free, free)],
                                        stiff_s[np.ix_(free, pres)] @ eps[pres])
            errors["free strain"] = np.max(np.abs(eps[free] - eps_free)) / eps_scale
        # 4. Prescribed strains against the expected engineering strains.
        if args.expect:
            expect = np.array(args.expect)
            if expect.size != 3:
                print("FAIL  --expect needs three values e11,e22,e33")
                return 1
            errors["target"] = (np.max(np.abs(eps[pres] - expect[pres]))
                                / np.max(np.abs(expect[pres])))

        # Tolerances of each check.
        limits = {"shear strain": args.rtol, "shear/free stress": args.stress_tol,
                  "hooke": args.rtol, "free strain": args.rtol, "target": args.target_rtol}
        for name, err in errors.items():
            worst[name] = max(worst[name], err)
            if err > limits[name]:
                n_fail += 1
                print(f"FAIL  elt {i_elt}: {name} error {err:.3g} (limit {limits[name]:g}); "
                      f"strain {eps}, stress {sig}")
                break

    summary = ", ".join(f"{name} {err:.3g}" for name, err in worst.items())
    print(f"check_triaxial_elastic: step {args.step}, {len(stress)} elements, {n_fail} failed; "
          f"worst: {summary}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
