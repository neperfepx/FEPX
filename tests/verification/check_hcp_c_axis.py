#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check an hcp single crystal under uniaxial stress along its c-axis (independent c33).

Source : FEPX verification tests (ROADMAP step 6.6: hcp with c33 independent of the other
         constants, D3, which FEPX v2 cannot represent). Development tool (CLAUDE.md §6).

The crystal's c-axis is the loading axis z (identity orientation), the lateral faces are free
(uniaxial_minimal), and the strengths are constant (hardening none). Three checks:

  1. elastic step (--elastic-step, a small strain --elastic-strain): the stress is uniaxial and
     sigma_zz = eps_zz / S33, S = C^-1 the compliance of the hexagonal stiffness (c11, c12, c13,
     c33, c44, c66 = (c11 - c12)/2; Voigt form, so S33 is that of the normal components), to
     --tol-elastic relative (eps_zz the true strain log(1 + engineering strain)); with v2's
     constraint c33 = c11 + c12 - c13 instead of the given c33 the check fails (titanium's
     c33 = 180.7 GPa against 185.4 GPa: a difference of a few percent in sigma_zz, far above
     the tolerance; both values are printed);
  2. the stress stays uniaxial at --flow-step (check_uniaxial_stress.py);
  3. flow stress at --flow-step: basal and prismatic systems have no resolved shear stress
     under c-axis loading, the 12 pyramidal <c+a> systems {10-11}<11-23> the same Schmid factor
         S = cos(phi) cos(lambda),  n = (1, 1/sqrt3, a/c) for (10-11), d = (3/2, 3 sqrt3/2, 3 c/a)
         for [11-23] (Miller-Bravais to Cartesian, a1 along x, c along z),
     so check_flow_stress.py applies with N = 12 (the steady flow stress depends on the slip
     strength, m and gammadot_0 only, not on the elastic constants).

Usage: check_hcp_c_axis.py SIM_DIR --c11 . --c12 . --c13 . --c33 . --c44 . --c-over-a .
                           --elastic-step N --elastic-strain E --flow-step N --flow-strain E
                           --crss C --m M --gammadot0 G0 --rate R [--tol-elastic T] [--tol-flow T]
Exit status 0 if every check passes, 1 otherwise.

Functions:
    compliance_33   S33 of the hexagonal stiffness
    schmid_pyramidal  Schmid factor of {10-11}<11-23> under c-axis loading
    main            command-line interface
"""

import argparse
import math
import os
import sys

import numpy as np

import check_flow_stress
import check_uniaxial_stress


def compliance_33(c11, c12, c13, c33, c44):
    """Return S33 of the hexagonal stiffness (Voigt 6x6, c66 = (c11 - c12)/2)."""
    stiff = np.zeros((6, 6))
    stiff[:3, :3] = [[c11, c12, c13], [c12, c11, c13], [c13, c13, c33]]
    stiff[3, 3] = stiff[4, 4] = c44
    stiff[5, 5] = 0.5 * (c11 - c12)
    return float(np.linalg.inv(stiff)[2, 2])


def schmid_pyramidal(c_over_a):
    """Return the Schmid factor of {10-11}<11-23> under loading along the c-axis."""
    normal = np.array([1.0, 1.0 / math.sqrt(3.0), 1.0 / c_over_a])
    direction = np.array([1.5, 1.5 * math.sqrt(3.0), 3.0 * c_over_a])
    return float(normal[2] / np.linalg.norm(normal) * direction[2] / np.linalg.norm(direction))


def main(argv=None):
    """Command-line interface; returns the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("sim_dir")
    for name in ("c11", "c12", "c13", "c33", "c44", "c-over-a", "elastic-strain", "flow-strain",
                 "crss", "m", "gammadot0", "rate"):
        parser.add_argument(f"--{name}", type=float, required=True)
    for name in ("elastic-step", "flow-step"):
        parser.add_argument(f"--{name}", type=int, required=True)
    parser.add_argument("--tol-elastic", type=float, default=1.0e-4)
    parser.add_argument("--tol-flow", type=float, default=1.0e-5)
    args = parser.parse_args(argv)
    ok = True

    # 1. Elastic step: uniaxial, sigma_zz = eps_zz / S33.
    ok &= check_uniaxial_stress.main([args.sim_dir, "--step", str(args.elastic_step),
                                      "--tol", "1e-6"]) == 0
    stress = np.loadtxt(os.path.join(args.sim_dir, "results", "elts", "stress",
                                     f"stress.step{args.elastic_step}"), ndmin=2)
    expected = math.log(1.0 + args.elastic_strain) / compliance_33(
        args.c11, args.c12, args.c13, args.c33, args.c44)
    constrained = math.log(1.0 + args.elastic_strain) / compliance_33(
        args.c11, args.c12, args.c13, args.c11 + args.c12 - args.c13, args.c44)
    worst = float(np.max(np.abs(stress[:, 2] - expected))) / abs(expected)
    good = worst <= args.tol_elastic
    print(f"check_hcp_c_axis: step {args.elastic_step}, sigma_zz = eps/S33 = {expected:.6e} "
          f"(with v2's constrained c33: {constrained:.6e}): largest relative difference "
          f"{worst:.3e} (tol {args.tol_elastic:.1e}) -> {'ok' if good else 'FAIL'}")
    ok &= good

    # 2-3. Uniaxial at the flow step, and the flow stress of the 12 pyramidal <c+a> systems.
    ok &= check_uniaxial_stress.main([args.sim_dir, "--step", str(args.flow_step),
                                      "--tol", "1e-6"]) == 0
    schmid = schmid_pyramidal(args.c_over_a)
    print(f"check_hcp_c_axis: Schmid factor of {{10-11}}<11-23> along c: {schmid:.15f}")
    ok &= check_flow_stress.main([args.sim_dir, "--step", str(args.flow_step),
                                  "--crss", str(args.crss),
                                  "--m", str(args.m), "--gammadot0", str(args.gammadot0),
                                  "--rate", str(args.rate), "--strain", str(args.flow_strain),
                                  "--schmid", repr(schmid), "--num-systems", "12",
                                  "--tol", str(args.tol_flow)]) == 0
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
