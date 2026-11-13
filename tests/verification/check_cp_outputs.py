#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check the crystal-plasticity outputs of a run for consistency with each other (6.5, 6.5b).

Source : FEPX verification tests (ROADMAP steps 6.5, 6.5b). Development tool (CLAUDE.md §6).

The outputs are written at full precision (17 digits) at the end of every step; with one
increment per step (from --first-step on; step 1 may contain the initialization increment),
the values at steps n - 1 and n are those of the start and end of an increment of length
dt = t_n - t_(n-1). For every element:

  1. slip:          slip_n - slip_(n-1) = dt sliprate_n (backward Euler, signed);
  2. strain_pl:     accumulated in the crystal frame, written in the sample frame:
                    R_n^T E_n R_n - R_(n-1)^T E_(n-1) R_(n-1) = dt R_n^T D_n R_n,
                    E = strain_pl, D = defrate_pl, R = rot_c2s of the written orientation;
  3. strain_pl_eq:  sqrt(2/3 E':E');
  4. sliprate:      the power law of the written rss and crss,
                    gammadot_0 |rss/crss|^(1/m) sgn(rss) (systems above the cutoff of
                    slip_kinetics_mod, |rss/crss| > 2e-6 here; below it the rate must be tiny);
  5. crss (with --h0, --tau-s, --tau-0; saturation_isotropic, n = 1): the backward-Euler
                    recursion crss_n = (crss_(n-1) + c tau_s)/(1 + c), c = dt h_0 sum|gammadot_n|/
                    (tau_s - tau_0), the same for every system;
  6. rotrate (with --rotrate): check_rotrate.py's finite difference of the orientations.
Step 6.5b (every step from --first-step; d = defrate, D^p = defrate_pl, s = stress):
  7. strain:        strain_n - strain_(n-1) = dt d_n (sample frame); strain_eq = sqrt(2/3 e':e');
  8. strain_pl_sample: its increment = dt D^p_n (sample frame); strain_pl_sample_eq likewise;
  9. velgrad:       sym(L) = d; defrate_eq = sqrt(2/3 d':d'); defrate_pl_eq = sqrt(2/3 D^p:D^p);
 10. stress_eq = sqrt(3/2 s':s'); strain_el_eq = sqrt(2/3 e_el':e_el');
 11. workrate = s : d, workrate_pl = s : D^p, work_pl_n - work_pl_(n-1) =
                    dt/2 (workrate_pl_n + workrate_pl_(n-1));
 12. rotrate_slip = (-w23, w13, -w12) of spinrate (w12, w13, w23), to --rtol-spin (default
                    1e-3: the strengths of the last two state iterations differ by less than
                    toler_state tau_0 = 1e-5 tau_0, so the slip rates by less than about
                    (1/m) 1e-5 = 2e-4 relative for m = 0.05): rotrate_slip is
                    the slip spin of the last orientation update (slip rates before the last
                    strength update, as v2's d_rstar_slip), spinrate that of the written slip
                    rates (recomputed with the final strengths, consistent with crss and rss);
 13. rotrate = rotrate_spin - rotrate_slip, to --rtol-rotrate: the two parts are those of the
                    orientation update (v2 d_rstar_spin, d_rstar_slip, d_rstar);
 14. convergence (results/convergence): v2's header, then rows "incr iter SA|NR six-norms
                    cg_iter"; every increment from 0 (initialization) to the last has rows
                    numbered 1, 2, ..., and the method is SA or NR.

The orientations must be Rodrigues vectors, passive convention (check_rotrate.py). Tolerances
(relative to the largest value of the compared quantity over the elements of the step):
--tol for checks 1-4 and 7-11 (default 1e-9), --rtol-rotrate for checks 6 and 13 (default 1e-9),
--rtol-spin for check 12 (default 1e-3); check 5 is
absolute, --tol-hard times tau_0 (default 1e-5: the backward-Euler strength update takes a
Newton step before its convergence test, as v2, which is exact for n = 1, so the strengths
satisfy the recursion up to
the state-loop tolerance 1e-5 of the slip rates; measured 1e-7 to 3.1e-6 tau_0).

Usage: check_cp_outputs.py SIM_DIR --times T1 T2 ... --m M --gammadot0 G0 [--first-step K]
                           [--h0 H --tau-s S --tau-0 T] [--rotrate] [--tol T] [--tol-hard T]
                           [--rtol-rotrate T] [--rtol-spin T]
Exit status 0 if every check passes, 1 otherwise.

Functions:
    sym             symmetric 3x3 tensors of rows 11 22 33 23 13 12
    report          print one check's result
    equiv           equivalent value sqrt(c a':a')
    scalar          one value per element of a scalar output
    check_derived   checks 7-13 at one step
    check_convergence  check 14
    main            command-line interface
"""

import argparse
import sys

import numpy as np

from check_rotrate import read_rows, rot_c2s_from_rodrigues_passive, rotation_vector


def sym(rows):
    """Return the symmetric tensors (n, 3, 3) of rows of tensor components 11 22 33 23 13 12."""
    mat = np.zeros((len(rows), 3, 3))
    mat[:, 0, 0], mat[:, 1, 1], mat[:, 2, 2] = rows[:, 0], rows[:, 1], rows[:, 2]
    mat[:, 1, 2] = mat[:, 2, 1] = rows[:, 3]
    mat[:, 0, 2] = mat[:, 2, 0] = rows[:, 4]
    mat[:, 0, 1] = mat[:, 1, 0] = rows[:, 5]
    return mat


def report(name, step, diff, scale, tol):
    """Print a check result (difference relative to scale against tol); return pass/fail."""
    rel = diff / scale if scale > 0.0 else diff
    good = rel <= tol
    print(f"check_cp_outputs: step {step}, {name}: largest difference {diff:.3e} "
          f"(relative {rel:.3e}; tol {tol:.1e}) -> {'ok' if good else 'FAIL'}")
    return good


def equiv(mat, factor):
    """Return sqrt(factor a':a') of symmetric tensors (n, 3, 3)."""
    dev = mat - np.einsum("eii->e", mat)[:, None, None] / 3.0 * np.eye(3)
    return np.sqrt(factor * np.einsum("eij,eij->e", dev, dev))


def scalar(sim, name, step):
    """Return the one value per element of a scalar output."""
    return read_rows(sim, name, step)[:, 0]


def check_derived(sim, step, dt, args):
    """Checks 7-13 of the module docstring at one step; return pass/fail."""
    ok = True
    defrate = sym(read_rows(sim, "defrate", step))
    defrate_pl = sym(read_rows(sim, "defrate_pl", step))
    stress = sym(read_rows(sim, "stress", step))

    # 7-8. Total strain and sample-frame plastic strain: sample-frame sums; equivalents.
    for name, rate in (("strain", defrate), ("strain_pl_sample", defrate_pl)):
        cur = sym(read_rows(sim, name, step))
        incr = cur - sym(read_rows(sim, name, step - 1))
        ok &= report(f"{name} increment = dt rate", step, np.max(np.abs(incr - dt * rate)),
                     np.max(np.abs(dt * rate)), args.tol)
        expected = equiv(cur, 2.0 / 3.0)
        ok &= report(f"{name}_eq = sqrt(2/3 e':e')", step,
                     np.max(np.abs(scalar(sim, name + "_eq", step) - expected)),
                     np.max(expected), args.tol)

    # 9. Velocity gradient and the equivalent rates.
    velgrad = read_rows(sim, "velgrad", step).reshape(-1, 3, 3)
    sym_part = 0.5 * (velgrad + np.transpose(velgrad, (0, 2, 1)))
    ok &= report("sym(velgrad) = defrate", step, np.max(np.abs(sym_part - defrate)),
                 np.max(np.abs(defrate)), args.tol)
    for name, rate in (("defrate_eq", defrate), ("defrate_pl_eq", defrate_pl)):
        expected = equiv(rate, 2.0 / 3.0)
        ok &= report(f"{name} = sqrt(2/3 d':d')", step,
                     np.max(np.abs(scalar(sim, name, step) - expected)), np.max(expected),
                     args.tol)

    # 10. Equivalent stress and elastic strain.
    expected = equiv(stress, 1.5)
    ok &= report("stress_eq = sqrt(3/2 s':s')", step,
                 np.max(np.abs(scalar(sim, "stress_eq", step) - expected)), np.max(expected),
                 args.tol)
    expected = equiv(sym(read_rows(sim, "strain_el", step)), 2.0 / 3.0)
    ok &= report("strain_el_eq = sqrt(2/3 e':e')", step,
                 np.max(np.abs(scalar(sim, "strain_el_eq", step) - expected)), np.max(expected),
                 args.tol)

    # 11. Work rates and the plastic work (trapezoidal rule).
    for name, rate in (("workrate", defrate), ("workrate_pl", defrate_pl)):
        expected = np.einsum("eij,eij->e", stress, rate)
        ok &= report(f"{name} = stress : rate", step,
                     np.max(np.abs(scalar(sim, name, step) - expected)),
                     np.max(np.abs(expected)), args.tol)
    incr = scalar(sim, "work_pl", step) - scalar(sim, "work_pl", step - 1)
    expected = 0.5 * dt * (scalar(sim, "workrate_pl", step)
                           + scalar(sim, "workrate_pl", step - 1))
    ok &= report("work_pl increment = trapezoidal rule", step, np.max(np.abs(incr - expected)),
                 np.max(np.abs(expected)), args.tol)

    # 12-13. Slip and spin parts of the rotation rate.
    spin = read_rows(sim, "spinrate", step)
    slip = read_rows(sim, "rotrate_slip", step)
    axial = np.stack([-spin[:, 2], spin[:, 1], -spin[:, 0]], axis=1)
    ok &= report("rotrate_slip = axial vector of spinrate", step, np.max(np.abs(slip - axial)),
                 np.max(np.abs(axial)), args.rtol_spin)
    total = read_rows(sim, "rotrate", step)
    parts = read_rows(sim, "rotrate_spin", step) - slip
    ok &= report("rotrate = rotrate_spin - rotrate_slip", step, np.max(np.abs(total - parts)),
                 np.max(np.abs(total)), args.rtol_rotrate)
    return ok


def check_convergence(sim, last_incr):
    """Check 14 of the module docstring; return pass/fail."""
    with open(f"{sim}/results/convergence", encoding="ascii") as handle:
        lines = handle.read().splitlines()
    good = lines[0].split()[:4] == ["%", "incr", "iter", "NR"]
    expected_iter = {}
    for line in lines[1:]:
        fields = line.split()
        incr, it = int(fields[0]), int(fields[1])
        good &= len(fields) == 10 and fields[2] in ("SA", "NR")
        good &= it == expected_iter.get(incr, 0) + 1
        expected_iter[incr] = it
    good &= sorted(expected_iter) == list(range(0, last_incr + 1))
    print(f"check_convergence: {len(lines) - 1} rows, increments "
          f"{min(expected_iter)}..{max(expected_iter)} -> {'ok' if good else 'FAIL'}")
    return good


def main(argv=None):
    """Command-line interface; returns the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("sim_dir")
    parser.add_argument("--times", type=float, nargs="+", required=True)
    parser.add_argument("--first-step", type=int, default=1)
    parser.add_argument("--m", type=float, required=True)
    parser.add_argument("--gammadot0", type=float, required=True)
    parser.add_argument("--h0", type=float)
    parser.add_argument("--tau-s", type=float)
    parser.add_argument("--tau-0", type=float)
    parser.add_argument("--rotrate", action="store_true")
    parser.add_argument("--tol", type=float, default=1.0e-9)
    parser.add_argument("--tol-hard", type=float, default=1.0e-5)
    parser.add_argument("--rtol-rotrate", type=float, default=1.0e-9)
    parser.add_argument("--rtol-spin", type=float, default=1.0e-3)
    args = parser.parse_args(argv)

    times = [0.0] + args.times
    sim = args.sim_dir
    ok = True
    slipping = False
    for step in range(args.first_step, len(times)):
        dt = times[step] - times[step - 1]
        rot = rot_c2s_from_rodrigues_passive(read_rows(sim, "ori", step))
        rot_prev = rot_c2s_from_rodrigues_passive(read_rows(sim, "ori", step - 1))

        # 1. Accumulated slip.
        slip = read_rows(sim, "slip", step)
        slip_prev = read_rows(sim, "slip", step - 1)
        rate = read_rows(sim, "sliprate", step)
        incr = slip - slip_prev
        ok &= report("slip increment = dt sliprate", step, np.max(np.abs(incr - dt * rate)),
                     np.max(np.abs(incr)), args.tol)
        slipping = slipping or np.max(np.abs(rate)) > 0.0

        # 2. Plastic strain accumulated in the crystal frame.
        strain = sym(read_rows(sim, "strain_pl", step))
        strain_prev = sym(read_rows(sim, "strain_pl", step - 1))
        defrate = sym(read_rows(sim, "defrate_pl", step))
        lhs = (np.einsum("eji,ejk,ekl->eil", rot, strain, rot)
               - np.einsum("eji,ejk,ekl->eil", rot_prev, strain_prev, rot_prev))
        rhs = dt * np.einsum("eji,ejk,ekl->eil", rot, defrate, rot)
        ok &= report("strain_pl crystal-frame increment = dt defrate_pl", step,
                     np.max(np.abs(lhs - rhs)), np.max(np.abs(rhs)), args.tol)

        # 3. Equivalent plastic strain.
        dev = strain - np.einsum("eii->e", strain)[:, None, None] / 3.0 * np.eye(3)
        expected = np.sqrt(2.0 / 3.0 * np.einsum("eij,eij->e", dev, dev))
        written = read_rows(sim, "strain_pl_eq", step)[:, 0]
        ok &= report("strain_pl_eq = sqrt(2/3 E':E')", step,
                     np.max(np.abs(written - expected)), np.max(np.abs(expected)), args.tol)

        # 4. Slip rates from the written resolved shear stresses and strengths.
        rss = read_rows(sim, "rss", step)
        crss = read_rows(sim, "crss", step)
        ratio = rss / crss
        law = args.gammadot0 * np.abs(ratio) ** (1.0 / args.m) * np.sign(ratio)
        above = np.abs(ratio) > 2.0e-6
        diff = np.max(np.abs(np.where(above, rate - law, 0.0)))
        ok &= report("sliprate = power law of rss/crss", step, diff, np.max(np.abs(law)),
                     args.tol)
        ok &= report("sliprate below the cutoff", step,
                     np.max(np.abs(np.where(above, 0.0, rate)), initial=0.0),
                     np.max(np.abs(law)), args.tol)

        # 5. Isotropic hardening recursion (n = 1).
        if args.h0 is not None:
            crss_prev = read_rows(sim, "crss", step - 1)
            c = dt * args.h0 * np.sum(np.abs(rate), axis=1) / (args.tau_s - args.tau_0)
            expected = (crss_prev + (c * args.tau_s)[:, None]) / (1.0 + c)[:, None]
            ok &= report("crss backward-Euler recursion (absolute, / tau_0)", step,
                         np.max(np.abs(crss - expected)), args.tau_0, args.tol_hard)

        # 6. Lattice rotation rate against the orientations.
        if args.rotrate:
            fd = rotation_vector(np.einsum("eij,ekj->eik", rot, rot_prev)) / dt
            written = read_rows(sim, "rotrate", step)
            ok &= report("rotrate = rotvec(R R_prev^T)/dt", step,
                         np.max(np.abs(written - fd)), np.max(np.abs(fd)), args.rtol_rotrate)

        # 7-13. Step 6.5b outputs.
        ok &= check_derived(sim, step, dt, args)

    # 14. Convergence record.
    ok &= check_convergence(sim, len(times) - 1)

    if not slipping:
        print("check_cp_outputs: no system slips in the checked steps -> FAIL")
        ok = False
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
