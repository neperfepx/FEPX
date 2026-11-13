#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Kocks-Mecking hardening against Barton et al. (2013), Fig. 3 (step 7.2, approximate).

Source : dev/theory/addendum_kocks_mecking.md §9 item 6 (ROADMAP D60). Development tool
(CLAUDE.md §6).

Reads the results of tests/verification/cases/fepx_bcc_km_barton (unirradiated bcc iron with
Barton's parameters in FEPX units, slightly off [100], uniaxial at 100 /s to 4 % strain in 8
printed steps) and

  checks (pass/fail):
    1. at step 0 every strength is tau_0 + b mu_0 sqrt(h_n rho_0) with mu_0 Hill's average of
       iron's elastic constants (computed here): the network term b mu_0 sqrt(h_n rho_0) is
       32.08 MPa (relative 1e-9);
    2. at every printed step and element, crss - tau_0 = b mu_0 sqrt(h_n km_rho) for every system,
       with the written km_rho (relative 1e-9 of the strength);
    3. the element mean of the network term crss - tau_0 never decreases;
  reports (no pass/fail; the agreement is approximate, see below): the element means of the
  network term, km_rho, k_2 and the plastic strain at each step, and the network term
  interpolated at the plastic strains 0.002, 0.012, 0.022 and 0.032, next to Barton's Fig. 3 read at the axial strains 0.01, 0.02,
  0.03 and 0.04 (his crystal is irradiated and yields near 0.008, so these are taken as his
  plastic strains): about 40, 74, 95 and 112 MPa (read from the figure, +-3 MPa).

Not exact by construction: Barton's mu_0 is 86 GPa (Hill: 81.8 GPa, so the network term is 5 %
lower at equal km_rho), his crystal is irradiated (the damage term changes which systems slip, hence the total
slip rate per axial strain), his orientation is only "slightly off-[100]" and he has 48 systems
({123} included; FEPX bcc: 24), and his elastic constants are not stated.

Usage: check_km_barton.py SIM_DIR
Exit status 0 if every check passes, 1 otherwise.
"""

import argparse
import math
import os
import sys

import numpy as np

# Barton's parameters (FEPX units) and iron's elastic constants (the case's).
C11, C12, C44 = 231.4e3, 134.7e3, 116.4e3
TAU_0, KM_B, KM_HN, RHO_0 = 90.0, 2.48e-7, 0.125, 2.0e7
NUM_STEPS = 8
BARTON = [(0.002, 40.0), (0.012, 74.0), (0.022, 95.0), (0.032, 112.0)]


def hill_shear_modulus():
    """Hill's (1952) average of the Voigt and Reuss shear moduli of the cubic stiffness [MPa]."""
    g_voigt = (C11 - C12 + 3.0 * C44) / 5.0
    g_reuss = 5.0 / (4.0 / (C11 - C12) + 3.0 / C44)
    return 0.5 * (g_voigt + g_reuss)


def read(sim, var, step):
    """Element values of a variable at a step (one row per element)."""
    return np.atleast_2d(np.loadtxt(os.path.join(sim, "results", "elts", var,
                                                 f"{var}.step{step}")).T).T


def main(argv=None):
    """Command-line interface; returns the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("sim")
    args = parser.parse_args(argv)

    k_factor = KM_B * hill_shear_modulus() * math.sqrt(KM_HN)
    ok = True

    # 1. The initial strengths.
    crss0 = read(args.sim, "crss", 0)
    expected = TAU_0 + k_factor * math.sqrt(RHO_0)
    err = np.max(np.abs(crss0 / expected - 1.0))
    good = err < 1e-9
    ok = ok and good
    print(f"{'ok' if good else 'FAIL'}: step 0 strengths tau_0 + b mu_0 sqrt(h_n rho_0) = "
          f"{expected:.4f} MPa (network term {expected - TAU_0:.4f}; "
          f"mu_0 Hill {hill_shear_modulus():.1f} MPa); relative error {err:.1e}")

    # 2-3. Strengths from the density at every step; the mean network term per step.
    rows = []
    worst = 0.0
    for step in range(NUM_STEPS + 1):
        crss = read(args.sim, "crss", step)
        rho = read(args.sim, "km_rho", step)[:, 0]
        network_term = crss - TAU_0
        worst = max(worst, np.max(np.abs(network_term - k_factor * np.sqrt(rho)[:, None]) / crss))
        eps_p = read(args.sim, "strain_pl_eq", step)[:, 0] if step > 0 else np.zeros(len(rho))
        k_2 = read(args.sim, "km_k2", step)[:, 0] if step > 0 else np.zeros(len(rho))
        rows.append((0.005 * step, eps_p.mean(), network_term.mean(), rho.mean(), k_2.mean()))
    good = worst < 1e-9
    ok = ok and good
    print(f"{'ok' if good else 'FAIL'}: crss - tau_0 = b mu_0 sqrt(h_n km_rho) at every step, element "
          f"and system (worst relative {worst:.1e})")
    means = [r[2] for r in rows]
    good = all(b >= a - 1e-12 for a, b in zip(means, means[1:]))
    ok = ok and good
    print(f"{'ok' if good else 'FAIL'}: the mean network term never decreases")

    # Report.
    print("report: strain  plastic strain  mean network term [MPa]  mean km_rho [1/mm^2]  "
          "mean k_2")
    for strain, eps_p, network_term, rho, k_2 in rows:
        print(f"report: {strain:6.3f}  {eps_p:14.5f}  {network_term:14.2f}  {rho:19.4e}  "
              f"{k_2:8.3f}")
    eps = [r[1] for r in rows]
    print("report: plastic strain  FEPX network term [MPa]  Barton Fig. 3 [MPa]  difference")
    for eps_target, barton in BARTON:
        network_term = float(np.interp(eps_target, eps, means))
        print(f"report: {eps_target:14.3f}  {network_term:15.1f}  {barton:19.0f}  "
              f"{100.0 * (network_term / barton - 1.0):+9.0f} %")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
