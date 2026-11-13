#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check triaxial load control (FEPX steps 8.2b, 8.2c) on an isotropic elastic cube: the loads of
every increment on their prescribed values, and Hooke's law at the printed steps.

Source : FEPX verification tests (ROADMAP steps 8.2b, 8.2c,
         dev/theory/addendum_triaxial_load_control.md §3.1, §3.4 and §5 V1, V2). Development tool
         (CLAUDE.md §6: Python 3 + NumPy for comparing outputs in tests).

Theory (addendum §3.1, §3.4, iteration (A)): the control component p moves at a constant strain
rate; every load-controlled component j (rows --rows) has its face velocity iterated so that its
value q_j (face force F_j [N], column F<j> of results/forces/<j>1, or F_j over that face's current
area [MPa]) is within its tolerance of the value to reach, in every increment:
  path steps (the control has a load or stress target, v2 CSR):
      q_j^id = q_j,start + lambda (q_j* - q_j,start), lambda = (q_p - q_p,start)/(q_p* - q_p,start);
  time steps (the control has a strain or time target):
      q_j^id = q_j,start + (q_j* - q_j,start)(t - t_start)/(t_end - t_start),
with the start values those at the end of the previous step (step 1: 0, Q-T4). A path step ends
within its tolerance of the control target (landing, addendum §3.6). Step 8.2c: a row given as a
rate r (spec "r<value>") has q_j^id = q_j,start + r (t - t_start) in every step; a row held
(spec "h", e.g. the control in a dwell) q_j^id = q_j,start. With a load or stress rate of the
control (iteration (B), mode "clr") the control row is checked like the others, in time. On an isotropic elastic cube
with normal face conditions the solution is homogeneous: every element's normal stresses equal
the face forces over the current face areas, there are no shear stresses, and the elastic strain
e gives the Kirchhoff stress tau = lambda tr(e) I + 2 mu e, sigma = tau / det(I + e) (stress and
strain_el: tensor components 11 22 33 23 13 12).

Checks:
  1. every increment of every step: |q_j - q_j^id| <= tol_j (1e-9 relative slack; the control row
     too in "clr" steps);
  2. path steps: the control value at the step's last increment within its tolerance of the
     target;
  3. each step in --hooke-steps: the element stresses equal F_i/A_i on the faces with a force
     row (all three, --faces) and Hooke's law, to --rtol of max |sigma| (at least --atol, for
     a step that ends at zero stress).

Usage:
  check_triaxial_load_control.py SIM_DIR --control 11 --modes path path
      --rows 11:load:-300,-600 22:load:0,0 33:load:300,600 [--load-tol T] [--stress-tol T]
      --lame LAMBDA MU [--hooke-steps 1 2] [--rtol R] [--atol A]
  (a step's spec in --rows: a target, r<rate> for a rate row, h for a held value; --modes path,
  time or clr per step)
Exit status 0 when every check passes, 1 otherwise (message on standard output).
"""

import argparse
import os
import sys

import numpy as np

from check_uniaxial_elastic import read_step


def read_forces(sim_dir, comp):
    """Force table of face <axis>1: step, incr, force along comp, area, time per row."""
    face = "xyz"[comp - 1] + "1"
    table = np.loadtxt(os.path.join(sim_dir, "results", "forces", face), comments="%", ndmin=2)
    return table[:, 0].astype(int), table[:, 1 + comp], table[:, 5], table[:, 6]


def main(argv=None):
    """Command-line interface; returns the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("sim_dir")
    parser.add_argument("--control", type=int, choices=(11, 22, 33), required=True)
    parser.add_argument("--modes", nargs="+", choices=("path", "time", "clr"), required=True,
                        help="path, time or clr (the control row checked in time), per step")
    parser.add_argument("--rows", nargs="+", required=True,
                        help="comp:kind:target1,target2,... (kind load [N] or stress [MPa]); "
                             "the control's row for path steps (its value in time steps unused)")
    parser.add_argument("--load-tol", type=float, default=0.1)
    parser.add_argument("--stress-tol", type=float, default=0.1)
    parser.add_argument("--lame", type=float, nargs=2, required=True, help="lambda mu [MPa]")
    parser.add_argument("--hooke-steps", type=int, nargs="*", default=[])
    parser.add_argument("--rtol", type=float, default=1.0e-6)
    parser.add_argument("--atol", type=float, default=0.0, help="absolute floor [MPa]")
    args = parser.parse_args(argv)

    # Values of every row per increment.
    rows = {}
    steps = time = None
    for spec in args.rows:
        comp, kind, targets = spec.split(":")
        comp = int(comp)
        steps, force, area, time = read_forces(args.sim_dir, comp // 11)
        values = force / area if kind == "stress" else force
        tol = args.stress_tol if kind == "stress" else args.load_tol
        rows[comp] = (kind, values, tol, targets.split(","))
    num_steps = len(args.modes)

    failures = []
    worst = {comp: 0.0 for comp in rows}
    for k in range(1, num_steps + 1):
        idx = np.nonzero(steps == k)[0]
        if idx.size == 0:
            failures.append(f"step {k}: no increments")
            continue
        start = idx[0] - 1
        mode = args.modes[k - 1]
        if mode == "path":
            # Progress of the control along its step.
            _, q_p, tol_p, targets_p = rows[args.control]
            target_p = float(targets_p[k - 1])
            lam = (q_p[idx] - q_p[start]) / (target_p - q_p[start])
            # 2. The control lands.
            if abs(q_p[idx[-1]] - target_p) > tol_p:
                failures.append(f"step {k}: control ends at {q_p[idx[-1]]:.6e}, more than {tol_p} "
                                f"from {target_p}")
        else:
            lam = (time[idx] - time[start]) / (time[idx[-1]] - time[start])
        elapsed = time[idx] - time[start]
        # 1. Every load-controlled row on its values (the control too in clr steps).
        for comp, (kind, values, tol, targets) in rows.items():
            if comp == args.control and mode != "clr":
                continue
            spec = targets[k - 1]
            if spec == "h":
                ideal = values[start] + 0.0 * elapsed
            elif spec.startswith("r"):
                ideal = values[start] + float(spec[1:]) * elapsed
            else:
                ideal = values[start] + lam * (float(spec) - values[start])
            error = np.max(np.abs(values[idx] - ideal))
            worst[comp] = max(worst[comp], error)
            if error > tol * (1 + 1e-9):
                failures.append(f"step {k}, component {comp}: {kind} off its values by "
                                f"{error:.3e} > {tol}")
        print(f"step {k} ({mode}): {idx.size} increments")
    for comp, error in worst.items():
        print(f"component {comp}: largest |q - q_id| {error:.3e}")

    # 3. Homogeneous Hooke's law at the printed steps.
    lam_e, mu_e = args.lame
    for k in args.hooke_steps:
        stress = read_step(args.sim_dir, "stress", k)
        strain = read_step(args.sim_dir, "strain_el", k)
        last = np.nonzero(steps == k)[0][-1]
        face_stress = []
        for comp in (1, 2, 3):
            st, force, area, _ = read_forces(args.sim_dir, comp)
            face_stress.append(force[last] / area[last])
        scale = max(np.max(np.abs(stress)), args.atol / args.rtol)
        for e, (sig, eps) in enumerate(zip(stress, strain)):
            e_mat = np.array([[eps[0], eps[5], eps[4]], [eps[5], eps[1], eps[3]],
                              [eps[4], eps[3], eps[2]]])
            tau = lam_e * np.trace(e_mat) * np.eye(3) + 2.0 * mu_e * e_mat
            sigma = tau / np.linalg.det(np.eye(3) + e_mat)
            expect = np.array([sigma[0, 0], sigma[1, 1], sigma[2, 2], sigma[1, 2], sigma[0, 2],
                               sigma[0, 1]])
            if np.max(np.abs(sig - expect)) > args.rtol * scale:
                failures.append(f"step {k}, element {e + 1}: stress {sig} differs from Hooke's "
                                f"law {expect}")
                break
            if np.max(np.abs(sig[:3] - face_stress)) > args.rtol * scale or \
                    np.max(np.abs(sig[3:])) > args.rtol * scale:
                failures.append(f"step {k}, element {e + 1}: stress {sig} is not the face "
                                f"stresses {face_stress} (homogeneous, no shear)")
                break
        print(f"step {k}: Hooke's law and face stresses checked on {len(stress)} elements")

    for message in failures:
        print("FAIL: " + message)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
