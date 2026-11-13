#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check the time integration of the strengths, the stress and the orientation of the point model.

Source : dev/theory/addendum_hardening.md §5 test 3 (step 6.4), extended by
         dev/theory/addendum_point_verification.md (step 7.4a). Development tool (CLAUDE.md §6).

The input is written by tests/material/drive_hardening_point.f90 (cases in its header): a single
crystal point under a constant deformation rate {d} and spin w, from a zero state, integrated by
FEPX's point model (staggered forward/backward-Euler strengths, orientation by the exponential
map) with dt = 0.8 ... 0.05 s. This script integrates the rate equations independently by the
classical Runge-Kutta method (RK4) with a small step, its own slip systems and stiffness.

Cases 1-6 (cube orientation, uniaxial, fcc; no lattice spin by symmetry): the rate equations
(Mandel 11 22 33 23 13 12)

    {de/dt} = {d} - sum_a gammadot^a {p^a},   {tau} = [C]{e},
    gammadot^a = gammadot_0 (|tau^a|/crss^a)^(1/m) sgn(tau^a),   tau^a = {p^a}.{tau},
    isotropic:   dcrss^a/dt = h_0 ((tau_s - crss^1)/(tau_s - tau_0))^n sum_b |gammadot^b|
                 (every a),
    anisotropic: dcrss^a/dt = h_0 ((tau_s - crss^a)/(tau_s - tau_0))^n sum_b H_ab |gammadot^b|,
    Kocks-Mecking (cases 4 and 5, step 7.2, addendum_kocks_mecking.md §9 item 4):
                 d(km_rho)/dt = Gammadot (k_1 sqrt(km_rho) - k_2 km_rho),
                 k_2 = k_20 (gammadot_k0/Gammadot)^(1/n_km),
                 crss^a = tau_0 + b mu_0 sqrt(h_n km_rho),
                 with mu_0 Hill's average shear modulus of the stiffness (computed here), from
                 rho_0 = 2e7 (case 4) and 1e4 /mm^2 (case 5); integrated in km_rho itself (FEPX
                 integrates sqrt(km_rho)),
    Kocks-Mecking with irradiation softening (case 6, step 7.3, addendum_irradiation.md §10
                 item 4): also dH/dt = -eta sum_a (N^a:H) N^a |gammadot^a|, N^a = n^a n^a,
                 crss^a = tau_0 + b mu_0 (sqrt(h_n km_rho) + sqrt(h_d <N^a:H>_+)), from rho_0 = 2e7
                 /mm^2 and H_0 = 3 d_l rho_l (1/4) sum_<111> (I - n n) (rho_l 3.61e12 /mm^3,
                 d = 100 b, h_d 0.675, eta 100).

Cases 7-18 (addendum_point_verification.md §1-§4): the general equations with lattice spin, in
3x3 sample-frame tensors E (elastic strain), D, W (imposed), R = rot_c2s from the Rodrigues
vector r (P5), [A, B] = AB - BA:

    tau^a = p^a_c : (C_c (R^T E R)),   gammadot^a by the power law,
    Dp = R (sum_a gammadot^a p^a_c) R^T,   Wslip = R (sum_a gammadot^a q^a_c) R^T,
    dE/dt - 1/2 [E, [E, dE/dt]] = D - Dp - [E, W] + [E, [E, Dp]]                (P2),
    Wp = W - [E, Dp] - 1/2 [E, dE/dt]  (P1),   Omega = Wp - Wslip  (P3),
    omega = (Omega_32, Omega_13, Omega_21),   dr/dt = 1/2 (omega + (omega.r) r + omega x r) (P4),

the isotropic or anisotropic law (with saturation evolution tau_s = tau_s0 (Gammadot/gdot_s0)^m',
the cyclic modifier with dA^a/dt = |gammadot^a| and its clamp, reset of A at the reversal, or
the precipitate strength in tau_0), or constant strengths (the Kocks-Mecking limits 15-18).

Checks, per case, at each output time t > 0:
  1. the error of the finest FEPX run (dt 0.05 s) is below the tolerance: cases 1-6, the
     Kirchhoff difference tau_33 - tau_11 and the sorted strengths (relative; --tol, --tol-irr
     for case 6, whose softening with eta = 100 makes the first-order error constant about 4
     times that of case 4); cases 7-18, the full Kirchhoff stress ||tau - tau_ref||/||tau_ref||,
     the sorted strengths, and the misorientation angle between R and R_ref divided by the
     reference's lattice rotation angle since the start (--tol-spin);
  2. at the last output time the error halves with dt (first order): the observed order
     log2(err(dt)/err(dt/2)) of the two finest pairs lies in [0.8, 1.2], for each measure whose
     finest error is above round-off (1e-10);
  3. the RK4 reference agrees with a run with twice its step to 1 % of the tolerance;
  4. cases 15-18: the strengths are exactly (1e-12) tau_0 (15) or tau_0 + b mu_0 sqrt(h_n rho_0)
     (16-18), and the offset from tau_0 falls by sqrt(100) = 10 from each of 16, 17 to the next.
  5. cases 7-18 except the cyclic case 11: the Richardson extrapolation 2 y(0.05) - y(0.1) of
     the two finest runs (first-order error removed) agrees with RK4 to --tol-extrap (stress,
     strengths) and --tol-extrap-ori (orientation). This resolves terms far below the
     first-order error, e.g. -1/2 [E, dE/dt] of (P1) (FEPX's -[E, E_0]/(2 dt)): without it
     the extrapolated shear-case errors are 1.4e-6 (stress) and 2.8e-5 (orientation).
At t = 0 the start orientation of cases 7-18 must equal (P5) of r_0 = (0.1, 0.2, 0.3) to 1e-14
(this pins the driver's orientation convention). Before any comparison the script checks its
own orientation integration (addendum_point_verification.md §2): (a) (P5) equals the axis-angle
matrix, (b) RK4 on (P4) with a constant omega equals exp(t [omega]x) R(r_0), (c) a reversed
sign of omega x r fails (b).

Usage: check_hardening_integration.py FILE [--tol T] [--tol-irr T] [--tol-spin T]
                                           [--tol-extrap T] [--tol-extrap-ori T]
                                           [--ref-dt DT] [--ref-dt-spin DT] [--cases N ...]
Exit status 0 if every check passes, 1 otherwise.
"""

import argparse
import itertools
import math
import sys

import numpy as np

# Material and loading of cases 1-6 (as in drive_hardening_point.f90).
C11, C12, C44 = 168400.0, 121400.0, 75400.0
M_RATE, GDOT0, TAU_0, H0, TAU_S, N_EXP = 0.02, 1.0, 100.0, 200.0, 330.0, 1.5
LATENT = {1: None, 2: 1.4, 3: 5.0}
# Kocks-Mecking (Barton et al. 2013 in FEPX units) and the start densities of cases 4-6.
KM_B, KM_HN, KM_K1, KM_K20, KM_GK0, KM_N = 2.48e-7, 0.125, 4.5e5, 14.0, 1.0e10, 100.0
KM_RHO0 = {4: 2.0e7, 5: 1.0e4, 6: 2.0e7}
# Irradiation softening (case 6): h_d, eta, loop density [1/mm^3] and diameter [mm].
IRR_HD, IRR_ETA, IRR_LOOPS, IRR_DIAM = 0.675, 100.0, 3.61e12, 100.0 * 2.48e-7
D = 1.0e-3 * np.array([-0.5, -0.5, 1.0, 0.0, 0.0, 0.0])
S2 = math.sqrt(2.0)
# Cases 7-18: start Rodrigues vector, simple-shear rate [1/s], Kocks-Mecking limit densities.
ROD_0 = np.array([0.1, 0.2, 0.3])
GDOT_SHEAR = math.sqrt(3.0) * 1.0e-3
KM_RHO_LIMIT = {16: 1.0e4, 17: 1.0e2, 18: 1.0}
ROUND_OFF = 1.0e-10
# Smallest lattice rotation for which the relative orientation error is checked [deg].
MIN_LATTICE_DEG = 0.1
# Cyclic case: RK4 is first order across the active-set switches, so a finer reference step.
REF_DT_FACTOR = {11: 0.2}


def mandel(mat):
    """Mandel 6-vector (11 22 33 23 13 12) of a symmetric 3x3 matrix."""
    return np.array([mat[0, 0], mat[1, 1], mat[2, 2], S2 * mat[1, 2], S2 * mat[0, 2],
                     S2 * mat[0, 1]])


def from_mandel(vec):
    """Symmetric 3x3 matrix of a Mandel 6-vector."""
    a, b, c = vec[3] / S2, vec[4] / S2, vec[5] / S2
    return np.array([[vec[0], c, b], [c, vec[1], a], [b, a, vec[2]]])


def fcc_schmid():
    """Mandel symmetric Schmid tensors of the 12 fcc {111}<110> systems (any order and sign)."""
    out = []
    for normal in ([1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1]):
        n = np.array(normal, float) / math.sqrt(3.0)
        for s in itertools.product([1, 0, -1], repeat=3):
            s = np.array(s, float)
            if np.count_nonzero(s) == 2 and abs(s @ n) < 1e-12 and s[np.nonzero(s)[0][0]] > 0:
                s = s / S2
                out.append(mandel(0.5 * (np.outer(s, n) + np.outer(n, s))))
    assert len(out) == 12
    return np.array(out)


def fcc_plane_tensors():
    """Mandel N^a = n^a n^a of the slip-plane normals of the 12 systems, in fcc_schmid's order."""
    out = []
    for normal in ([1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1]):
        n = np.array(normal, float) / math.sqrt(3.0)
        out.extend([mandel(np.outer(n, n))] * 3)
    return np.array(out)


def irr_dmg_0():
    """H_0 (I4) of the four <111> loop normals, Mandel [1/mm^2]."""
    avg = np.zeros((3, 3))
    for normal in ([1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1]):
        n = np.array(normal, float) / math.sqrt(3.0)
        avg += (np.eye(3) - np.outer(n, n)) / 4.0
    return mandel(IRR_LOOPS * 3.0 * IRR_DIAM * avg)


def stiffness():
    """Cubic Mandel stiffness of copper [MPa]."""
    return cubic_stiffness(C11, C12, C44)


def cubic_stiffness(c11, c12, c44):
    """Cubic Mandel stiffness [MPa] (shear diagonal 2 C44)."""
    c = np.zeros((6, 6))
    c[:3, :3] = c12
    np.fill_diagonal(c[:3, :3], c11)
    c[3:, 3:] = 2.0 * c44 * np.eye(3)
    return c


def hexagonal_stiffness(c11, c12, c13, c33, c44):
    """Hexagonal Mandel stiffness, c-axis = crystal 3 [MPa] (shear 2C44, 2C44, 2C66)."""
    c = np.array([[c11, c12, c13], [c12, c11, c13], [c13, c13, c33]], float)
    out = np.zeros((6, 6))
    out[:3, :3] = c
    out[3:, 3:] = np.diag([2.0 * c44, 2.0 * c44, c11 - c12])
    return out


def hill_shear_modulus():
    """Hill's (1952) average of the Voigt and Reuss shear moduli of the cubic stiffness [MPa]."""
    g_voigt = (C11 - C12 + 3.0 * C44) / 5.0
    g_reuss = 5.0 / (4.0 / (C11 - C12) + 3.0 / C44)
    return 0.5 * (g_voigt + g_reuss)


def km_strengths(rho, dmg=None):
    """Kocks-Mecking strengths of the 12 systems from the density, with the defect term of H
    when given [MPa]."""
    crss = np.full(12, TAU_0 + KM_B * hill_shear_modulus() * math.sqrt(KM_HN * rho))
    if dmg is not None:
        contact = np.maximum(fcc_plane_tensors() @ dmg, 0.0)
        crss = crss + KM_B * hill_shear_modulus() * np.sqrt(IRR_HD * contact)
    return crss


def rates(y, p, c, law):
    """Rates of the state y = (e(6), crss(12)), or (e(6), km_rho) for Kocks-Mecking, or
    (e(6), km_rho, H(6)) with irradiation (cases 1-6)."""
    if law >= 4:
        e, rho = y[:6], y[6]
        dmg = y[7:13] if law == 6 else None
        crss = km_strengths(rho, dmg)
        rss = p @ (c @ e)
        gdot = GDOT0 * (np.abs(rss) / crss) ** (1.0 / M_RATE) * np.sign(rss)
        shrate = np.sum(np.abs(gdot))
        k2 = KM_K20 * (KM_GK0 / shrate) ** (1.0 / KM_N) if shrate > 0.0 else 0.0
        drho = shrate * (KM_K1 * math.sqrt(rho) - k2 * rho)
        if law == 6:
            planes = fcc_plane_tensors()
            ddmg = -IRR_ETA * ((planes @ dmg) * np.abs(gdot)) @ planes
            return np.concatenate([D - gdot @ p, [drho], ddmg])
        return np.concatenate([D - gdot @ p, [drho]])
    e, crss = y[:6], y[6:]
    tau = c @ e
    rss = p @ tau
    gdot = GDOT0 * (np.abs(rss) / crss) ** (1.0 / M_RATE) * np.sign(rss)
    de = D - gdot @ p
    if law == 1:
        ratio = max(TAU_S - crss[0], 0.0) / (TAU_S - TAU_0)
        dcrss = np.full(12, H0 * ratio ** N_EXP * np.sum(np.abs(gdot)))
    else:
        h = np.full((12, 12), LATENT[law])
        np.fill_diagonal(h, 1.0)
        ratio = np.maximum(TAU_S - crss, 0.0) / (TAU_S - TAU_0)
        dcrss = H0 * ratio ** N_EXP * (h @ np.abs(gdot))
    return np.concatenate([de, dcrss])


def rk4(law, times, step, p, c):
    """RK4 states at the given times (each a multiple of step; cases 1-6)."""
    y = np.concatenate([np.zeros(6), [KM_RHO0[law]] if law >= 4 else np.full(12, TAU_0)])
    if law == 6:
        y = np.concatenate([y, irr_dmg_0()])
    t, out = 0.0, {}
    for t_out in times:
        num = int(round((t_out - t) / step))
        for _ in range(num):
            k1 = rates(y, p, c, law)
            k2 = rates(y + 0.5 * step * k1, p, c, law)
            k3 = rates(y + 0.5 * step * k2, p, c, law)
            k4 = rates(y + step * k3, p, c, law)
            y = y + step / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        t = t_out
        out[t_out] = y.copy()
    return out


def measures(e, crss, c):
    """Kirchhoff tau_33 - tau_11 [MPa] and the sorted strengths [MPa]."""
    tau = c @ e
    return tau[2] - tau[0], np.sort(crss)


# ----------------------------------------------------------------------------------------------
# Cases 7-18 (addendum_point_verification.md).


def cross_matrix(v):
    """[v]x, the matrix of the cross product: [v]x u = v x u."""
    return np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])


def rot_from_rodrigues(r):
    """Active rotation of the Rodrigues vector r = n tan(phi/2) (P5, arXiv Eq. 21)."""
    rr = r @ r
    return ((1.0 - rr) * np.eye(3) + 2.0 * np.outer(r, r) + 2.0 * cross_matrix(r)) / (1.0 + rr)


def rot_from_axis_angle(vec):
    """exp([vec]x): the rotation by |vec| about vec/|vec| (Rodrigues' rotation formula)."""
    angle = np.linalg.norm(vec)
    if angle == 0.0:
        return np.eye(3)
    k = cross_matrix(vec / angle)
    return np.eye(3) + math.sin(angle) * k + (1.0 - math.cos(angle)) * (k @ k)


def rotation_angle(rot):
    """Angle of a rotation [rad], by atan2 (accurate for small angles)."""
    sin_part = 0.5 * np.linalg.norm([rot[2, 1] - rot[1, 2], rot[0, 2] - rot[2, 0],
                                     rot[1, 0] - rot[0, 1]])
    return math.atan2(sin_part, 0.5 * (np.trace(rot) - 1.0))


def rodrigues_rate(omega, r, sign=1.0):
    """dr/dt of (P4) (arXiv Eq. 31) for the sample-frame axial spin vector omega; sign = -1
    reverses omega x r (the mutation of self-check (c))."""
    return 0.5 * (omega + (omega @ r) * r + sign * np.cross(omega, r))


def self_check_orientation():
    """Self-checks (a)-(c) of addendum_point_verification.md §2; returns True if all pass."""
    rng = np.random.default_rng(7)
    ok_a = True
    for _ in range(20):
        r = rng.normal(size=3)
        phi = 2.0 * math.atan(np.linalg.norm(r))
        ok_a = ok_a and np.max(np.abs(rot_from_rodrigues(r) - rot_from_axis_angle(
            phi * r / np.linalg.norm(r)))) < 1e-14
    omega = np.array([0.3, -0.5, 0.7])
    errs = []
    for sign in (1.0, -1.0):
        r, h = ROD_0.copy(), 1.0e-3
        for _ in range(1000):
            k1 = rodrigues_rate(omega, r, sign)
            k2 = rodrigues_rate(omega, r + 0.5 * h * k1, sign)
            k3 = rodrigues_rate(omega, r + 0.5 * h * k2, sign)
            k4 = rodrigues_rate(omega, r + h * k3, sign)
            r = r + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        exact = rot_from_axis_angle(omega) @ rot_from_rodrigues(ROD_0)
        errs.append(np.max(np.abs(rot_from_rodrigues(r) - exact)))
    good = ok_a and errs[0] < 1e-12 and errs[1] > 1e-3
    print(f"{'ok' if good else 'FAIL'}: orientation self-checks: (a) (P5) = axis-angle "
          f"{'yes' if ok_a else 'no'}; (b) RK4 on (P4) vs exact {errs[0]:.1e}; "
          f"(c) reversed omega x r {errs[1]:.1e} (must fail)")
    return good


MANDEL_BASIS = np.array([from_mandel(row) for row in np.eye(6)])


def systems_from_pairs(pairs):
    """Mandel {p^a} (N, 6), 3x3 p^a and q^a (N, 3, 3) of (slip direction, plane normal) pairs."""
    p6, p3, q3 = [], [], []
    for s, n in pairs:
        s, n = s / np.linalg.norm(s), n / np.linalg.norm(n)
        assert abs(s @ n) < 1e-12
        sm = np.outer(s, n)
        p3.append(0.5 * (sm + sm.T))
        q3.append(0.5 * (sm - sm.T))
        p6.append(mandel(p3[-1]))
    return np.array(p6), np.array(p3), np.array(q3)


def unique_up_to_sign(vectors):
    """The vectors with one of each pair v, -v kept."""
    out = []
    for v in vectors:
        if not any(np.allclose(v, u) or np.allclose(v, -u) for u in out):
            out.append(v)
    return out


def cubic_pairs(plane_family, dir_family):
    """(direction, normal) pairs of a cubic {plane}<dir> family, from the permutations and signs
    of the Miller indices."""
    def family(idx):
        vecs = set()
        for perm in itertools.permutations(idx):
            for signs in itertools.product([1, -1], repeat=3):
                vecs.add(tuple(a * b for a, b in zip(perm, signs)))
        return unique_up_to_sign([np.array(v, float) for v in vecs])
    return [(s, n) for n in family(plane_family) for s in family(dir_family) if abs(s @ n) < 1e-12]


def hcp_pairs(c_over_a):
    """(direction, normal, family) of the hcp basal, prismatic and pyramidal <c+a> systems, from
    Miller-Bravais indices with the basis a1 = (1, 0, 0), a2 = (-1/2, sqrt3/2, 0),
    a3 = (-1/2, -sqrt3/2, 0), c = (0, 0, c/a): direction [u v t w] -> u a1 + v a2 + t a3 + w c;
    plane (h k i l) -> the normal h a1* + k a2* + l c* of the reciprocal basis of (a1, a2, c)."""
    a1, a2 = np.array([1.0, 0.0, 0.0]), np.array([-0.5, math.sqrt(3.0) / 2.0, 0.0])
    a3, cc = np.array([-0.5, -math.sqrt(3.0) / 2.0, 0.0]), np.array([0.0, 0.0, c_over_a])
    vol = a1 @ np.cross(a2, cc)
    a1s, a2s, cs = np.cross(a2, cc) / vol, np.cross(cc, a1) / vol, np.cross(a1, a2) / vol

    def direction(u, v, t, w):
        return u * a1 + v * a2 + t * a3 + w * cc

    def normal(h, k, i, l):
        return h * a1s + k * a2s + l * cs

    def mb_family(base):
        """Miller-Bravais indices from the permutations of the first three and signs."""
        out = set()
        for perm in itertools.permutations(base[:3]):
            for sgn in (1, -1):
                for lsgn in ((1, -1) if base[3] != 0 else (1,)):
                    out.add(tuple(sgn * x for x in perm) + (lsgn * base[3],))
        return sorted(out)

    out = []
    families = [((0, 0, 0, 1), (1, 1, -2, 0)), ((1, 0, -1, 0), (1, 1, -2, 0)),
                ((1, 0, -1, 1), (1, 1, -2, 3))]
    for fam, (plane, dirn) in enumerate(families):
        if plane == (0, 0, 0, 1):
            normals = [normal(0, 0, 0, 1)]
        else:
            normals = unique_up_to_sign([normal(*x) for x in mb_family(plane)])
        dirs = unique_up_to_sign([direction(*x) for x in mb_family(dirn)])
        for n in normals:
            for s in dirs:
                if abs(s @ n) < 1e-12:
                    out.append((s, n, fam))
    return out


class Case:
    """Material, law and loading of one of cases 7-18 (drive_hardening_point.f90 header)."""

    def __init__(self, num):
        self.num = num
        # Crystal and stiffness.
        if num == 13:
            pairs = cubic_pairs((1, 1, 0), (1, 1, 1)) + cubic_pairs((1, 1, 2), (1, 1, 1))
            fam = np.zeros(len(pairs), int)
            self.c_c = cubic_stiffness(231400.0, 134700.0, 116400.0)
        elif num == 14:
            trip = hcp_pairs(1.587)
            pairs, fam = [(s, n) for s, n, _ in trip], np.array([f for _, _, f in trip])
            self.c_c = hexagonal_stiffness(162400.0, 92000.0, 69000.0, 180700.0, 46700.0)
        else:
            pairs = cubic_pairs((1, 1, 1), (1, 1, 0))
            fam = np.zeros(len(pairs), int)
            self.c_c = cubic_stiffness(C11, C12, C44)
        self.p6, self.p3, self.q3 = systems_from_pairs(pairs)
        self.n = len(pairs)
        assert self.n == {13: 24, 14: 18}.get(num, 12)
        # Slip parameters per system; hcp per family; the precipitate strength in tau_0.
        self.m = np.full(self.n, M_RATE)
        self.tau_0 = np.full(self.n, TAU_0)
        if num == 14:
            self.tau_0 = np.array([150.0, 100.0, 250.0])[fam]
            self.m = np.array([0.02, 0.02, 0.05])[fam]
            # The isotropic law follows the first system (FEPX's: basal); put basal first.
            assert fam[0] == 0
        if num == 12:
            self.tau_0 = self.tau_0 + 100.0 * math.sqrt(0.01 * 1.0e-6 / 2.5e-7)
        # The law: 'iso', 'aniso' (latent 1.4) or 'const' (strengths fixed, cases 15-18).
        self.law = {9: "aniso"}.get(num, "const" if num >= 15 else "iso")
        self.sat_evol = num == 10
        self.cyclic = num == 11
        self.t_rev = 20.0 if num == 11 else None
        if num == 15:
            self.crss_const = self.tau_0.copy()
        elif num >= 16:
            self.crss_const = self.tau_0 + KM_B * hill_shear_modulus() * math.sqrt(
                KM_HN * KM_RHO_LIMIT[num])
        # Loading.
        if num in (8, 9):
            self.d = np.zeros((3, 3))
            self.d[0, 1] = self.d[1, 0] = 0.5 * GDOT_SHEAR
            self.w = np.zeros((3, 3))
            self.w[0, 1], self.w[1, 0] = 0.5 * GDOT_SHEAR, -0.5 * GDOT_SHEAR
        else:
            self.d = from_mandel(D)
            self.w = np.zeros((3, 3))
        dev = self.d - np.trace(self.d) / 3.0 * np.eye(3)
        self.d_eq = math.sqrt(2.0 / 3.0 * np.sum(dev * dev))

    def start(self):
        """Start state y = (e(6), r(3), crss(n) unless constant, A(n) if cyclic)."""
        parts = [np.zeros(6), ROD_0]
        if self.law != "const":
            parts.append(self.tau_0.copy())
        if self.cyclic:
            parts.append(np.zeros(self.n))
        return np.concatenate(parts)

    def split(self, y):
        """(E, r, crss, A) of a state."""
        n = self.n
        crss = self.crss_const if self.law == "const" else y[9:9 + n]
        acc = y[9 + n:9 + 2 * n] if self.cyclic else None
        return from_mandel(y[:6]), y[6:9], crss, acc

    def stress(self, y):
        """Sample-frame Kirchhoff stress (3x3) and R of a state."""
        e, r, _, _ = self.split(y)
        rot = rot_from_rodrigues(r)
        tau_c = from_mandel(self.c_c @ mandel(rot.T @ e @ rot))
        return rot @ tau_c @ rot.T, rot

    def rates(self, y, sign):
        """Rates of the state (addendum_point_verification.md §1-§3); sign = -1 after the
        reversal."""
        e, r, crss, acc = self.split(y)
        rot = rot_from_rodrigues(r)
        # Resolved shear stresses (frame-invariant: crystal frame) and slip rates.
        tau_c = self.c_c @ mandel(rot.T @ e @ rot)
        rss = self.p6 @ tau_c
        gdot = GDOT0 * (np.abs(rss) / crss) ** (1.0 / self.m) * np.sign(rss)
        # Plastic deformation rate and slip spin, sample frame.
        dp = rot @ np.tensordot(gdot, self.p3, 1) @ rot.T
        wslip = rot @ np.tensordot(gdot, self.q3, 1) @ rot.T
        d, w = sign * self.d, sign * self.w

        def comm(a, b):
            return a @ b - b @ a

        # (P2): the elastic strain rate, a 6x6 linear system built column by column.
        rhs = d - dp - comm(e, w) + comm(e, comm(e, dp))
        inner = e @ MANDEL_BASIS - MANDEL_BASIS @ e
        outer = MANDEL_BASIS - 0.5 * (e @ inner - inner @ e)
        lmat = np.stack([mandel(x) for x in outer], axis=1)
        edot = from_mandel(np.linalg.solve(lmat, mandel(rhs)))
        # (P1), (P3), (P4): the lattice spin and the Rodrigues rate.
        wp = w - comm(e, dp) - 0.5 * comm(e, edot)
        omega_mat = wp - wslip
        omega = np.array([omega_mat[2, 1], omega_mat[0, 2], omega_mat[1, 0]])
        parts = [mandel(edot), rodrigues_rate(omega, r)]
        # Hardening.
        if self.law != "const":
            shrate = np.sum(np.abs(gdot))
            tau_s = np.full(self.n, TAU_S)
            if self.sat_evol and shrate > 0.0:
                tau_s[:] = TAU_S * (shrate / 1.0e-2) ** 0.1
            if self.law == "iso":
                if self.cyclic:
                    crit = 5.0e-3 * abs(crss[0] / tau_s[0]) ** 2.0
                    crit = 0.0 if crit < 2.22e-16 else crit
                    drive = np.sum(np.abs(gdot)[acc >= crit])
                    drive = min(max(drive, 1.0e-6 * self.d_eq), 10.0 * self.d_eq)
                else:
                    drive = shrate
                ratio = max(tau_s[0] - crss[0], 0.0) / (tau_s[0] - self.tau_0[0])
                parts.append(np.full(self.n, H0 * ratio ** N_EXP * drive))
            else:
                h = np.full((self.n, self.n), 1.4)
                np.fill_diagonal(h, 1.0)
                ratio = np.maximum(tau_s - crss, 0.0) / (tau_s - self.tau_0)
                parts.append(H0 * ratio ** N_EXP * (h @ np.abs(gdot)))
        if self.cyclic:
            parts.append(np.abs(gdot))
        return np.concatenate(parts)

    def rk4(self, times, step):
        """RK4 states at the given times (multiples of step), with the reversal (sign change,
        accumulated slip reset) at t_rev."""
        y, t, out, sign = self.start(), 0.0, {0.0: self.start()}, 1.0
        stops = sorted(set(times) | ({self.t_rev} if self.t_rev else set()))
        for t_stop in stops:
            for _ in range(int(round((t_stop - t) / step))):
                k1 = self.rates(y, sign)
                k2 = self.rates(y + 0.5 * step * k1, sign)
                k3 = self.rates(y + 0.5 * step * k2, sign)
                k4 = self.rates(y + step * k3, sign)
                y = y + step / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            t = t_stop
            out[t] = y.copy()
            if self.t_rev and abs(t - self.t_rev) < 1e-9:
                sign = -1.0
                y[9 + self.n:] = 0.0
        return out

    def errors(self, y_ref, e6, rot, crss):
        """Relative errors of a FEPX state (e Mandel, R, crss) against an RK4 state: full
        Kirchhoff stress, sorted strengths, misorientation over the lattice rotation angle."""
        tau_ref, rot_ref = self.stress(y_ref)
        # FEPX's stress with its own orientation.
        tau = rot @ from_mandel(self.c_c @ mandel(rot.T @ from_mandel(e6) @ rot)) @ rot.T
        err_s = np.linalg.norm(tau - tau_ref) / np.linalg.norm(tau_ref)
        crss_ref = self.split(y_ref)[2]
        err_crss = np.max(np.abs(np.sort(crss) / np.sort(crss_ref) - 1.0))
        err_r = rotation_angle(rot_ref.T @ rot) / self.lattice_angle(y_ref)
        return err_s, err_crss, err_r

    def lattice_angle(self, y_ref):
        """Lattice rotation angle of a state since the start [rad]."""
        return rotation_angle(rot_from_rodrigues(y_ref[6:9]) @ rot_from_rodrigues(ROD_0).T)


def check_general(rows, num, tol, ref_dt, args_richardson):
    """Checks 1-5 of one of cases 7-18; returns True if all pass. args_richardson: the
    tolerances of check 5 (stress and strengths, orientation)."""
    case = Case(num)
    mine = rows[rows[:, 0] == num]
    times = sorted(t for t in set(mine[:, 2]) if t > 0.0)
    dts = sorted(set(mine[:, 1]), reverse=True)
    ok = True
    # Start orientation (pins the driver's convention).
    start = mine[mine[:, 2] == 0.0][0]
    start_err = np.max(np.abs(start[11:20].reshape(3, 3, order="F") - rot_from_rodrigues(ROD_0)))
    if start_err > 1e-14:
        print(f"FAIL: case {num}: start orientation differs from (P5) of r_0 by {start_err:.1e}")
        ok = False
    ref_dt = ref_dt * REF_DT_FACTOR.get(num, 1.0)
    ref = case.rk4(times, ref_dt)
    ref2 = case.rk4(times, 2.0 * ref_dt)
    for t in times:
        e_r2 = case.errors(ref[t], ref2[t][:6], case.stress(ref2[t])[1], case.split(ref2[t])[2])
        ref_err = max(e_r2)
        errs, its = [], 0
        for dt in dts:
            row = mine[(mine[:, 1] == dt) & (mine[:, 2] == t)][0]
            n = int(row[4])
            errs.append(case.errors(ref[t], row[5:11], row[11:20].reshape(3, 3, order="F"),
                                    row[20:20 + n]))
            its = max(its, int(row[3]))
        errs = np.array(errs)
        # The orientation error only where the lattice has rotated enough to measure it.
        checked = [True, True, math.degrees(case.lattice_angle(ref[t])) > MIN_LATTICE_DEG]
        good = all(errs[-1, j] <= tol for j in range(3) if checked[j]) and ref_err < 0.01 * tol
        orders = np.full((len(dts) - 1, 3), np.nan)
        for j in range(3):
            if errs[-1, j] > ROUND_OFF:
                orders[:, j] = np.log2(errs[:-1, j] / errs[1:, j])
                if t == times[-1] and checked[j]:
                    if case.cyclic:
                        # Switches of the active set at increment boundaries make the order
                        # irregular: the mean order over the whole range instead.
                        mean = math.log2(errs[0, j] / errs[-1, j]) / (len(dts) - 1)
                        good = good and 0.8 <= mean <= 1.2
                    else:
                        good = good and all(0.8 <= o <= 1.2 for o in orders[-2:, j])
        # Check 5: the Richardson extrapolation of the two finest runs (not cyclic).
        extra = ""
        if t == times[-1] and not case.cyclic:
            rows_2 = [mine[(mine[:, 1] == dt) & (mine[:, 2] == t)][0] for dt in dts[-2:]]
            e_x = richardson_errors(case, ref[t], *rows_2)
            good = good and max(e_x[:2]) <= args_richardson[0] and e_x[2] <= args_richardson[1]
            extra = (f"; extrapolated (2 y(dt/2) - y(dt)) stress {e_x[0]:.1e}, strengths "
                     f"{e_x[1]:.1e}, orientation {e_x[2]:.1e}")
        # Cases 15-18: the strengths are exactly the constants.
        if case.law == "const":
            last = mine[(mine[:, 1] == dts[-1]) & (mine[:, 2] == t)][0]
            crss_err = np.max(np.abs(last[20:20 + case.n] / case.crss_const - 1.0))
            good = good and crss_err < 1e-12
        ok = ok and good

        def fmt(col):
            return ", ".join(f"{x:.2e}" for x in errs[:, col])

        def fmt_o(col):
            return ", ".join("-" if math.isnan(o) else f"{o:.2f}" for o in orders[:, col])
        print(f"{'ok' if good else 'FAIL'}: case {num} t {t:g} s: relative errors "
              f"(dt {', '.join(f'{d:g}' for d in dts)}): stress {fmt(0)} (orders {fmt_o(0)}); "
              f"strengths {fmt(1)} (orders {fmt_o(1)}); orientation {fmt(2)} "
              f"(orders {fmt_o(2)}; lattice rotation "
              f"{math.degrees(case.lattice_angle(ref[t])):.2f} deg"
              f"{'' if checked[2] else ', not checked'}); "
              f"reference {ref_err:.1e}; state iterations <= {its}{extra}")
    return ok


def richardson_errors(case, y_ref, row_coarse, row_fine):
    """Errors (as Case.errors) of the Richardson extrapolation 2 y(dt/2) - y(dt) of two FEPX
    runs: elastic strain and strengths linearly, the orientation along the geodesic,
    R = R_f exp(log(R_c^T R_f))."""
    n = case.n
    rot_c = row_coarse[11:20].reshape(3, 3, order="F")
    rot_f = row_fine[11:20].reshape(3, 3, order="F")
    step = rot_c.T @ rot_f
    angle = rotation_angle(step)
    axis = np.array([step[2, 1] - step[1, 2], step[0, 2] - step[2, 0], step[1, 0] - step[0, 1]])
    if angle > 0.0:
        axis = axis / np.linalg.norm(axis)
    rot = rot_f @ rot_from_axis_angle(angle * axis)
    e6 = 2.0 * row_fine[5:11] - row_coarse[5:11]
    crss = 2.0 * row_fine[20:20 + n] - row_coarse[20:20 + n]
    return case.errors(y_ref, e6, rot, crss)


def check_km_limit_offsets(rows):
    """Check 4 (cases 16-18): the offsets of the constant strengths from tau_0 scale as
    sqrt(rho_0), i.e. fall by 10 from case to case."""
    offsets = []
    for num in (16, 17, 18):
        row = rows[(rows[:, 0] == num) & (rows[:, 2] == 40.0)][-1]
        offsets.append(row[20] - TAU_0)
    ratios = [offsets[0] / offsets[1], offsets[1] / offsets[2]]
    good = all(abs(q / 10.0 - 1.0) < 1e-10 for q in ratios)
    print(f"{'ok' if good else 'FAIL'}: Kocks-Mecking limit: offsets from tau_0 "
          f"{', '.join(f'{o:.6e}' for o in offsets)} MPa (rho_0 1e4, 1e2, 1 /mm^2); ratios "
          f"{', '.join(f'{q:.12f}' for q in ratios)} (sqrt(100) = 10)")
    return good


def check_cube(rows, law, args, p, c):
    """Checks 1-3 of one of cases 1-6 (the original check)."""
    mine = rows[(rows[:, 0] == law) & (rows[:, 2] > 0.0)]
    times = sorted(set(mine[:, 2]))
    dts = sorted(set(mine[:, 1]), reverse=True)
    ref = rk4(law, times, args.ref_dt, p, c)
    ref2 = rk4(law, times, 2.0 * args.ref_dt, p, c)
    ok = True
    for t in times:
        if law == 6:
            crss_of = lambda y: km_strengths(y[6], y[7:13])
        elif law >= 4:
            crss_of = lambda y: km_strengths(y[6])
        else:
            crss_of = lambda y: y[6:]
        s_ref, crss_ref = measures(ref[t][:6], crss_of(ref[t]), c)
        s_ref2, crss_ref2 = measures(ref2[t][:6], crss_of(ref2[t]), c)
        ref_err = max(abs(s_ref2 / s_ref - 1.0), np.max(np.abs(crss_ref2 / crss_ref - 1.0)))
        errs_s, errs_crss, its = [], [], 0
        for dt in dts:
            row = mine[(mine[:, 1] == dt) & (mine[:, 2] == t)][0]
            s, crss = measures(row[5:11], row[20:32], c)
            errs_s.append(abs(s / s_ref - 1.0))
            errs_crss.append(np.max(np.abs(crss / crss_ref - 1.0)))
            its = max(its, int(row[3]))
        orders = [math.log2(errs_crss[i] / errs_crss[i + 1]) for i in range(len(dts) - 1)]
        tol = args.tol_irr if law == 6 else args.tol
        good = max(errs_s[-1], errs_crss[-1]) <= tol and ref_err < 0.01 * tol
        if t == times[-1]:
            good = good and all(0.8 <= o <= 1.2 for o in orders[-2:])
        ok = ok and good
        print(f"{'ok' if good else 'FAIL'}: case {law} t {t:g} s: relative errors "
              f"(dt {', '.join(f'{d:g}' for d in dts)}): "
              f"tau33-tau11 {', '.join(f'{x:.2e}' for x in errs_s)}; "
              f"strengths {', '.join(f'{x:.2e}' for x in errs_crss)}; "
              f"orders {', '.join(f'{o:.2f}' for o in orders)}; "
              f"reference {ref_err:.1e}; state iterations <= {its}")
    return ok


def main(argv=None):
    """Command-line interface; returns the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("file")
    parser.add_argument("--tol", type=float, default=2.0e-4)
    parser.add_argument("--tol-irr", type=float, default=1.0e-3)
    parser.add_argument("--tol-spin", type=float, default=1.0e-3)
    parser.add_argument("--tol-extrap", type=float, default=5.0e-7)
    parser.add_argument("--tol-extrap-ori", type=float, default=5.0e-6)
    parser.add_argument("--ref-dt", type=float, default=2.0e-3)
    parser.add_argument("--ref-dt-spin", type=float, default=1.0e-2)
    parser.add_argument("--cases", type=int, nargs="*", default=list(range(1, 19)))
    args = parser.parse_args(argv)

    rows = np.atleast_2d(np.loadtxt(args.file))
    p, c = fcc_schmid(), stiffness()
    ok = self_check_orientation()
    for num in args.cases:
        if num <= 6:
            ok = check_cube(rows, num, args, p, c) and ok
        else:
            ok = check_general(rows, num, args.tol_spin, args.ref_dt_spin,
                               (args.tol_extrap, args.tol_extrap_ori)) and ok
    if all(num in args.cases for num in (16, 17, 18)):
        ok = check_km_limit_offsets(rows) and ok
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
