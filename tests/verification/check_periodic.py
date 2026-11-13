#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check an FEPX periodic run (ROADMAP step 8.3b; dev/theory/addendum_periodic_mpc.md §5).

Source : FEPX verification tests (step 8.3b, V1-V3, V8). Development tool (CLAUDE.md §6).

Reads results/macro (step incr time, Fbar 11 12 13 21 ... 33, Fbar-dot (same order), stress 11
22 33 23 13 12), the node coordinates and velocities, the element stresses (and velgrad) of the
listed steps, and the mesh's $Periodicity section (inputs/simulation.msh). Checks at each step N
(the last increment of the step):
  relations   every relation (s, p, P): v_s - v_p = Fbar-dot (P o l) and x_s - x_p = Fbar (P o l)
              within --rtol of the largest jump (l from the mesh's drive pairs: the first relation
              with vector +-e_j, as FEPX);
  --homogeneous TOL  every element stress within TOL of the field scale of the mean (a uniform
              field), and the mean equal to the volume-averaged stress of results/macro (V8);
  --imposed ij=value ...  Fbar-dot_ij equal to the value (relative 1e-10 of the largest);
  --symmetric  Fbar-dot symmetric (small strain);
  --free c ...  the volume-averaged stress components (11 22 33 23 13 12) at most --stol of the
              largest one (free macroscopic components are stress free: exact for the Cauchy
              stress only while no shear stress acts, see --free-nominal);
  --free-nominal kj ...  the nominal stress P = det(Fbar) sigma Fbar^-T of the volume-averaged
              stress (results/macro) at most --stol of its largest component, for the free
              components kj of Fbar-dot (11 12 13 21 ... 33): P is the stress conjugate to
              Fbar-dot, which a free component's equation sets to zero (addendum_periodic_mpc.md
              §2: the sum of the drive pair's relation reactions); with --symmetric (small
              strain: one unknown per shear pair) a shear kj is checked as P_kj + P_jk. Under
              shear the Cauchy components of free normal directions are then of the order of
              the shear stress times the strain, not zero (step 9.6);
  --elastic C11 C12  isotropic elasticity (C44 = (C11 - C12)/2): the volume-averaged stress
              equals C : e with e = sym(Fbar - I), within --erel of the largest stress (small
              strains);
  --velgrad   every element's velocity gradient equals Fbar-dot Fbar^-1 within --rtol (the
              homogeneous solution of V1 with the nine components);
  --poisson NU  the lateral rates of L = Fbar-dot Fbar^-1: L22 = L33 = -NU L11 (relative 1e-4).

Usage: check_periodic.py SIM_DIR --steps N ... [options above] [--rtol R]
Exit status 0 if every check passes, 1 otherwise.
"""

import argparse
import os
import re
import sys

import numpy as np

STRESS_COMPONENTS = ["11", "22", "33", "23", "13", "12"]
FORTRAN_E3 = re.compile(r"^([+-]?(?:\d+\.?\d*|\.\d+))([+-]\d{3})$")


def number(token):
    """A Fortran real, also with a three-digit exponent written without the E."""
    match = FORTRAN_E3.match(token)
    if match:
        return float(match.group(1) + "e" + match.group(2))
    return float(token)


def read_table(path):
    """Rows of numbers of a result file (comment lines starting with % skipped)."""
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.strip() and not line.startswith("%"):
                rows.append([number(t) for t in line.split()])
    return np.array(rows)


def read_mesh(path):
    """Initial coordinates (array by id - 1) and periodicity relations (s, p, P)."""
    with open(path) as handle:
        text = handle.read()
    nodes = text.split("$Nodes\n", 1)[1].split("$EndNodes", 1)[0].split("\n")[1:]
    coo = np.array([[float(v) for v in line.split()[1:4]] for line in nodes if line.strip()])
    lines = text.split("$Periodicity\n", 1)[1].split("$EndPeriodicity", 1)[0].split("\n")[1:]
    relations = np.array([[int(v) for v in line.split()] for line in lines if line.strip()])
    return coo, relations


def lengths_of(coo, relations):
    """Domain lengths from the drive pairs (first relation with vector +-e_j)."""
    lengths = np.zeros(3)
    for s, p, *vector in relations:
        vector = np.array(vector)
        if np.abs(vector).sum() == 1:
            j = int(np.argmax(np.abs(vector)))
            if lengths[j] == 0.0:
                lengths[j] = abs(coo[s - 1, j] - coo[p - 1, j])
    return lengths


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sim_dir")
    parser.add_argument("--steps", type=int, nargs="+", required=True)
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    parser.add_argument("--homogeneous", type=float)
    parser.add_argument("--imposed", nargs="*", default=[])
    parser.add_argument("--symmetric", action="store_true")
    parser.add_argument("--free", nargs="*", default=[], choices=STRESS_COMPONENTS)
    parser.add_argument("--free-nominal", nargs="*", default=[],
                        choices=[f"{k}{j}" for k in "123" for j in "123"])
    parser.add_argument("--stol", type=float, default=1.0e-8)
    parser.add_argument("--elastic", type=float, nargs=2)
    parser.add_argument("--erel", type=float, default=1.0e-4)
    parser.add_argument("--velgrad", action="store_true")
    parser.add_argument("--poisson", type=float)
    args = parser.parse_args()

    results = os.path.join(args.sim_dir, "results")
    coo0, relations = read_mesh(os.path.join(args.sim_dir, "inputs", "simulation.msh"))
    lengths = lengths_of(coo0, relations)
    macro = read_table(os.path.join(results, "macro"))
    failed = False

    def report(ok, text):
        nonlocal failed
        print(("ok    " if ok else "FAIL  ") + text)
        failed = failed or not ok

    for step in args.steps:
        row = macro[macro[:, 0] == step][-1]
        fbar = row[3:12].reshape(3, 3)
        fdot = row[12:21].reshape(3, 3)
        stress = row[21:27]
        jumps = (relations[:, 2:] * lengths)
        sec = relations[:, 0] - 1
        pri = relations[:, 1] - 1

        # Relations on the velocities and the coordinates.
        vel = read_table(os.path.join(results, "nodes", "vel", f"vel.step{step}"))
        coo = read_table(os.path.join(results, "nodes", "coo", f"coo.step{step}"))
        for name, field, grad in (("vel", vel, fdot), ("coo", coo, fbar)):
            actual = field[sec] - field[pri]
            misfit = np.abs(actual - jumps @ grad.T).max() / np.abs(actual).max()
            report(misfit <= args.rtol, f"step {step}: relations ({name}) misfit {misfit:.3e}")

        # Imposed components and symmetry.
        for item in args.imposed:
            comp, value = item.split("=")
            i, j = int(comp[0]) - 1, int(comp[1]) - 1
            err = abs(fdot[i, j] - float(value)) / np.abs(fdot).max()
            report(err <= 1e-10, f"step {step}: Fbar-dot {comp} = {fdot[i, j]:.6e} (given {value})")
        if args.symmetric:
            err = np.abs(fdot - fdot.T).max() / np.abs(fdot).max()
            report(err <= 1e-10, f"step {step}: Fbar-dot symmetric ({err:.3e})")

        # Homogeneous field and the volume average.
        if args.homogeneous is not None:
            elt = read_table(os.path.join(results, "elts", "stress", f"stress.step{step}"))
            mean = elt.mean(axis=0)
            spread = np.abs(elt - mean).max() / np.abs(mean).max()
            report(spread <= args.homogeneous, f"step {step}: homogeneous stress ({spread:.3e})")
            err = np.abs(mean - stress).max() / np.abs(mean).max()
            report(err <= args.homogeneous, f"step {step}: element stresses = macro stress "
                   f"({err:.3e})")

        # Free components.
        scale = np.abs(stress).max()
        for comp in args.free:
            value = stress[STRESS_COMPONENTS.index(comp)]
            report(abs(value) <= args.stol * scale,
                   f"step {step}: free stress {comp} = {value:.3e} (scale {scale:.3e})")

        # Free components of Fbar-dot: the conjugate nominal stress P = det(F) sigma F^-T.
        if args.free_nominal:
            sigma = np.array([[stress[0], stress[5], stress[4]],
                              [stress[5], stress[1], stress[3]],
                              [stress[4], stress[3], stress[2]]])
            nominal = np.linalg.det(fbar) * sigma @ np.linalg.inv(fbar).T
            nominal_scale = np.abs(nominal).max()
            for comp in args.free_nominal:
                k, j = int(comp[0]) - 1, int(comp[1]) - 1
                value = nominal[k, j]
                # Small strain: the shear pair is one unknown, conjugate to P_kj + P_jk.
                if args.symmetric and k != j:
                    value = nominal[k, j] + nominal[j, k]
                report(abs(value) <= args.stol * nominal_scale,
                       f"step {step}: free nominal stress {comp} = {value:.3e} "
                       f"(scale {nominal_scale:.3e})")

        # Isotropic elasticity at small strain.
        if args.elastic:
            c11, c12 = args.elastic
            c44 = 0.5 * (c11 - c12)
            eps = 0.5 * (fbar - np.eye(3) + (fbar - np.eye(3)).T)
            sig = c12 * np.trace(eps) * np.eye(3) + 2.0 * c44 * eps
            expected = np.array([sig[0, 0], sig[1, 1], sig[2, 2], sig[1, 2], sig[0, 2],
                                 sig[0, 1]])
            err = np.abs(stress - expected).max() / np.abs(expected).max()
            report(err <= args.erel, f"step {step}: stress = C : e ({err:.3e})")

        # Velocity gradient of the homogeneous solution, and the lateral rates.
        lbar = fdot @ np.linalg.inv(fbar)
        if args.velgrad:
            grads = read_table(os.path.join(results, "elts", "velgrad", f"velgrad.step{step}"))
            err = np.abs(grads - lbar.reshape(9)).max() / np.abs(lbar).max()
            report(err <= args.rtol, f"step {step}: element velgrad = Fbar-dot Fbar^-1 "
                   f"({err:.3e})")
        if args.poisson is not None:
            for k in (1, 2):
                ratio = -lbar[k, k] / lbar[0, 0]
                report(abs(ratio - args.poisson) <= 1e-4 * args.poisson,
                       f"step {step}: -L{k + 1}{k + 1}/L11 = {ratio:.6f} (nu {args.poisson})")

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
