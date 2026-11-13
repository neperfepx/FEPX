#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check the Kocks-Mecking and irradiation outputs of a run against each other and the laws (9.6).

Source : FEPX verification tests (ROADMAP step 9.6, D88). Development tool (CLAUDE.md §6).
Theory : addendum_kocks_mecking.md (K2), (K2s); addendum_irradiation.md (I3), §4 (integration).

The outputs are written at full precision at the end of every step; with one increment per step
(from --first-step on), the values at steps n - 1 and n are those of the start and end of an
increment of length dt = t_n - t_(n-1). For every element (Gammadot = sum |sliprate_n|):

  1. km_k2:       k_2 = k_20 (gammadot_k0 / Gammadot)^(1/n_k), zero at Gammadot = 0;
  2. km_rhorate:  rho dot = Gammadot (k_1 sqrt(rho_n) - k_2 rho_n)   (K2) at the written values;
  3. km_rho:      backward Euler in s = sqrt(rho) (addendum_kocks_mecking.md (K2s)):
                  sqrt(rho_n) - sqrt(rho_(n-1)) = dt rhodot_n / (2 sqrt(rho_n));
  4. dmg:         backward Euler (addendum_irradiation.md §4): H_n - H_(n-1) = dt Hdot_n.

Checks 1 and 2 relate values written from the same state: tolerance --tol (default 1e-12)
relative to the largest value of the step. Checks 3 and 4 relate the update made in the point
solve to the rates written at the slip rates of the final stress, which differ by the state-loop
tolerance: tolerance --tol-update (default 1e-4) relative to the largest increment of the step.
At least one increment must have Gammadot > 0 and a non-zero H increment (the checks would
otherwise be empty).

Usage: check_km_irr_outputs.py SIM_DIR --times T1 T2 ... --k1 K1 --k20 K20 --gammadot-k0 G
                               --n-k N [--first-step K] [--tol T] [--tol-update T]
Exit status 0 if every check passes, 1 otherwise.

Functions:
    report  print one check's result
    main    command-line interface
"""

import argparse
import sys

import numpy as np

from check_rotrate import read_rows


def report(name, step, diff, scale, tol):
    """Print one check: name, step, max |difference|, its ratio to scale; return pass/fail.

    name   (str)   what is checked
    step   (int)   step number
    diff   (float) largest absolute difference over the elements
    scale  (float) the reference magnitude (largest value or increment of the step)
    tol    (float) allowed diff / scale
    ok     (bool)  diff <= tol scale
    """
    ratio = diff / scale if scale > 0.0 else diff
    ok = ratio <= tol
    print(f"{'ok  ' if ok else 'FAIL'} step {step}: {name}: max |diff| {diff:.3e}, "
          f"relative {ratio:.3e} (tolerance {tol:.1e})")
    return ok


def main(argv=None):
    """Command-line interface; returns the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("sim_dir")
    parser.add_argument("--times", type=float, nargs="+", required=True)
    parser.add_argument("--first-step", type=int, default=1)
    parser.add_argument("--k1", type=float, required=True)
    parser.add_argument("--k20", type=float, required=True)
    parser.add_argument("--gammadot-k0", type=float, required=True)
    parser.add_argument("--n-k", type=float, required=True)
    parser.add_argument("--tol", type=float, default=1.0e-12)
    parser.add_argument("--tol-update", type=float, default=1.0e-4)
    args = parser.parse_args(argv)

    times = [0.0] + args.times
    sim = args.sim_dir
    ok = True
    plastic = False
    for step in range(args.first_step, len(times)):
        dt = times[step] - times[step - 1]
        shrate = np.sum(np.abs(read_rows(sim, "sliprate", step)), axis=1)
        rho = read_rows(sim, "km_rho", step)[:, 0]
        rho_prev = read_rows(sim, "km_rho", step - 1)[:, 0]
        k2 = read_rows(sim, "km_k2", step)[:, 0]
        rhorate = read_rows(sim, "km_rhorate", step)[:, 0]
        dmg = read_rows(sim, "dmg", step)
        dmg_prev = read_rows(sim, "dmg", step - 1)
        dmgrate = read_rows(sim, "dmgrate", step)

        # 1. Recovery coefficient k_2 of the written slip rates.
        k2_law = np.zeros_like(shrate)
        active = shrate > 0.0
        k2_law[active] = args.k20 * (args.gammadot_k0 / shrate[active]) ** (1.0 / args.n_k)
        ok &= report("km_k2 = k20 (gammadot_k0/Gammadot)^(1/n)", step,
                     np.max(np.abs(k2 - k2_law)), np.max(np.abs(k2_law)), args.tol)

        # 2. Density rate (K2) at the written density, k_2 and slip rates.
        rate_law = shrate * (args.k1 * np.sqrt(rho) - k2 * rho)
        ok &= report("km_rhorate = Gammadot (k1 sqrt(rho) - k2 rho)", step,
                     np.max(np.abs(rhorate - rate_law)), np.max(np.abs(rate_law)), args.tol)

        # 3. Backward Euler of s = sqrt(rho).
        s_incr = np.sqrt(rho) - np.sqrt(rho_prev)
        s_rhs = dt * rhorate / (2.0 * np.sqrt(rho))
        ok &= report("sqrt(rho) increment = dt km_rhorate / (2 sqrt(rho))", step,
                     np.max(np.abs(s_incr - s_rhs)), np.max(np.abs(s_incr)), args.tol_update)

        # 4. Backward Euler of the damage tensor.
        h_incr = dmg - dmg_prev
        ok &= report("dmg increment = dt dmgrate", step,
                     np.max(np.abs(h_incr - dt * dmgrate)), np.max(np.abs(h_incr)),
                     args.tol_update)
        plastic = plastic or (np.max(shrate) > 0.0 and np.max(np.abs(h_incr)) > 0.0)

    # The checks are not empty: some increment slipped and changed the damage.
    if not plastic:
        print("FAIL: no increment with slip and a damage increment (the checks are empty)")
        ok = False
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
