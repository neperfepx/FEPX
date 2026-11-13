# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Write the yield-surface vertex fixtures of tests/material/test_yield_vertices.f90 with KD25's
own implementation (Kasemer & Dawson 2025, github.com/mkasemer/scys; not part of FEPX).

Usage: python3 -I kd25_dump_vertices.py <scys clone> > kd25_vertices.txt   (README.md)

Each case: a header line "<structure> <c_over_a> <num_slip> <num_vertices>", a line with the
strength of each slip system (KD25's system order, which is FEPX's for the cases below), and one
line per vertex: the deviatoric 5-vector of KD25 Eq. 8 (FEPX addendum_5vectors.md (F1)).
"""
import contextlib
import io
import sys

import numpy as np

sys.path.insert(0, sys.argv[1])
import scys_calc_scys_vertices as kd25  # noqa: E402

# KD25's bct system 26 has the normal (1 2 1) with the direction [1 0 1], which are not
# orthogonal (n . d = 2): not a slip system. FEPX and v2 have (-1 2 1)[1 0 1] there
# (crystal_mod bct table). The 5-vector of the corrected system replaces it (README.md).
_schmid_vecs = kd25.scys_calc_schmid_vecs


def _schmid_vecs_corrected(ctype, covera=1.633):
    pdev = _schmid_vecs(ctype, covera)
    if ctype == "bct":
        n = np.array([-1.0, 2.0, 1.0]) / np.sqrt(6.0)
        d = np.array([1.0, 0.0, 1.0]) / np.sqrt(2.0)
        p = 0.5 * (np.outer(d, n) + np.outer(n, d))
        pdev[:, 25] = [np.sqrt(0.5) * (p[0, 0] - p[1, 1]), np.sqrt(1.5) * p[2, 2],
                       np.sqrt(2.0) * p[1, 2], np.sqrt(2.0) * p[0, 2], np.sqrt(2.0) * p[0, 1]]
    return pdev


kd25.scys_calc_schmid_vecs = _schmid_vecs_corrected
scys_calc_scys_vertices = kd25.scys_calc_scys_vertices

CASES = [
    # fcc {111}<110>, unit strengths (KD25 §5.1: 56 vertices).
    ("fcc", 1.0, [1.0]*12),
    # bcc {110}<111> with system 1 5 % weaker (KD25 §6.2: 116 vertices).
    ("bcc", 1.0, [0.95] + [1.0]*11),
    # hcp, c/a 1.587, basal 1, prismatic 1.25, pyramidal <c+a> 1.7 (KD25's example).
    ("hcp", 1.587, [1.0]*3 + [1.25]*3 + [1.7]*12),
    # hcp, c/a 1.59, pyramidal <c+a> only (basal and prismatic out of reach; KD25 §5.2: 92).
    ("hcp", 1.59, [100.0]*6 + [1.0]*12),
    # bct, unit strengths (KD25's bct ignores c/a: c/a = 1).
    ("bct", 1.0, [1.0]*32),
]

for structure, c_over_a, strengths in CASES:
    tau = np.array([[t, t] for t in strengths])
    with contextlib.redirect_stdout(io.StringIO()):
        verts, num_vertices, _ = scys_calc_scys_vertices(structure, tau, covera=c_over_a)
    print(f"{structure} {c_over_a:.6f} {len(strengths)} {num_vertices}")
    print(" ".join(f"{t:.6f}" for t in strengths))
    for v in verts.T:
        print(" ".join(f"{x:.6f}" for x in v))
