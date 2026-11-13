#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Unit tests of tests/tools/compare_sim.py (ROADMAP step 0.3).

Each test builds small synthetic .sim trees in a temporary directory and checks that the
comparison passes or fails as specified in the compare_sim.py docstring and
tests/regression/README.md. Run directly or through CTest (test tools_compare_sim).
"""

import os
import sys
import tempfile
import unittest

# Import compare_sim from this directory (tests/tools/).
TOOLS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(TOOLS_DIR))
import compare_sim  # noqa: E402

# Reference stress field: 3 elements x 6 components, with near-zero shear components.
REF_STRESS = [
    [100.0, 50.0, 1000.0, 1.0e-3, -2.0e-3, 0.0],
    [110.0, 55.0, 1010.0, 3.0e-3, 1.0e-3, 0.0],
    [-90.0, 45.0, 990.0, -1.0e-3, 0.0, 0.0],
]


def write_sim(root, files):
    """Create root/results/<path> for each {path: rows} entry; rows are lists of floats."""
    for rel, rows in files.items():
        path = os.path.join(root, "results", rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="ascii") as handle:
            for row in rows:
                handle.write(" ".join(f"{v:.7E}" for v in row) + "\n")


def scaled(rows, factor):
    """Return rows with every value multiplied by factor."""
    return [[v * factor for v in row] for row in rows]


class CompareSimTest(unittest.TestCase):
    """Pass/fail behaviour of compare_sim on synthetic result trees."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.ref = os.path.join(self.tmp.name, "ref.sim")
        self.test = os.path.join(self.tmp.name, "test.sim")

    def tearDown(self):
        self.tmp.cleanup()

    def run_compare(self, options=None):
        """Return (number of failed files, results)."""
        results = compare_sim.compare_sim(self.test, self.ref, options or compare_sim.Options())
        return sum(not r.passed for r in results), results

    def test_identical_trees_pass(self):
        files = {"elts/stress/stress.step1": REF_STRESS, "nodes/coo/coo.step1": [[0, 0, 1]]}
        write_sim(self.ref, files)
        write_sim(self.test, files)
        n_fail, results = self.run_compare()
        self.assertEqual(n_fail, 0)
        self.assertEqual(len(results), 2)

    def test_relative_difference_within_rtol_passes(self):
        write_sim(self.ref, {"elts/stress/stress.step1": REF_STRESS})
        write_sim(self.test, {"elts/stress/stress.step1": scaled(REF_STRESS, 1.0 + 5.0e-6)})
        self.assertEqual(self.run_compare()[0], 0)

    def test_relative_difference_beyond_rtol_fails(self):
        write_sim(self.ref, {"elts/stress/stress.step1": REF_STRESS})
        bad = [row[:] for row in REF_STRESS]
        bad[1][0] *= 1.0 + 1.0e-3
        write_sim(self.test, {"elts/stress/stress.step1": bad})
        n_fail, results = self.run_compare()
        self.assertEqual(n_fail, 1)
        self.assertIn("line 2, column 1", results[0].message)

    def test_every_column_is_compared(self):
        # v2's own --diff compared only the first value of each line (CHANGES.md C5).
        write_sim(self.ref, {"elts/stress/stress.step1": REF_STRESS})
        bad = [row[:] for row in REF_STRESS]
        bad[2][2] += 1.0
        write_sim(self.test, {"elts/stress/stress.step1": bad})
        n_fail, results = self.run_compare()
        self.assertEqual(n_fail, 1)
        self.assertIn("line 3, column 3", results[0].message)

    def test_small_component_uses_field_scale_floor(self):
        # A shear component of 1e-3 changed by 5e-4 (50 % relative) is within
        # ftol * scale = 1e-6 * 1010 ~ 1.0e-3, so it passes.
        write_sim(self.ref, {"elts/stress/stress.step1": REF_STRESS})
        near = [row[:] for row in REF_STRESS]
        near[0][3] += 5.0e-4
        write_sim(self.test, {"elts/stress/stress.step1": near})
        self.assertEqual(self.run_compare()[0], 0)

    def test_small_component_beyond_floor_fails(self):
        write_sim(self.ref, {"elts/stress/stress.step1": REF_STRESS})
        far = [row[:] for row in REF_STRESS]
        far[0][3] += 5.0e-3
        write_sim(self.test, {"elts/stress/stress.step1": far})
        self.assertEqual(self.run_compare()[0], 1)

    def test_exact_zero_field_requires_zero(self):
        write_sim(self.ref, {"nodes/disp/disp.step0": [[0.0, 0.0, 0.0]]})
        write_sim(self.test, {"nodes/disp/disp.step0": [[0.0, 1.0e-30, 0.0]]})
        self.assertEqual(self.run_compare()[0], 1)
        # ... unless an absolute tolerance is given.
        opts = compare_sim.Options(compare_sim.Tolerance(atol=1.0e-20))
        self.assertEqual(self.run_compare(opts)[0], 0)

    def test_missing_and_extra_files_fail(self):
        write_sim(self.ref, {"elts/stress/stress.step1": REF_STRESS,
                             "elts/stress/stress.step2": REF_STRESS})
        write_sim(self.test, {"elts/stress/stress.step1": REF_STRESS,
                              "elts/strain/strain.step1": REF_STRESS})
        n_fail, results = self.run_compare()
        self.assertEqual(n_fail, 2)
        messages = {r.path: r.message for r in results if not r.passed}
        self.assertIn("missing", messages["elts/stress/stress.step2"])
        self.assertIn("extra", messages["elts/strain/strain.step1"])

    def test_ignore_patterns(self):
        write_sim(self.ref, {"elts/stress/stress.step1": REF_STRESS,
                             "mesh/stress/stress.step1": REF_STRESS})
        write_sim(self.test, {"elts/stress/stress.step1": REF_STRESS})
        opts = compare_sim.Options(ignore=["mesh/*"])
        self.assertEqual(self.run_compare(opts)[0], 0)

    def test_shape_mismatch_fails(self):
        write_sim(self.ref, {"elts/stress/stress.step1": REF_STRESS})
        write_sim(self.test, {"elts/stress/stress.step1": REF_STRESS[:2]})
        n_fail, results = self.run_compare()
        self.assertEqual(n_fail, 1)
        self.assertIn("line count", results[0].message)
        write_sim(self.test, {"elts/stress/stress.step1": [r[:5] for r in REF_STRESS]})
        n_fail, results = self.run_compare()
        self.assertEqual(n_fail, 1)
        self.assertIn("value count", results[0].message)

    def test_nan_fails(self):
        write_sim(self.ref, {"elts/stress/stress.step1": REF_STRESS})
        bad = [row[:] for row in REF_STRESS]
        bad[0][0] = float("nan")
        write_sim(self.test, {"elts/stress/stress.step1": bad})
        n_fail, results = self.run_compare()
        self.assertEqual(n_fail, 1)
        self.assertIn("non-finite", results[0].message)

    def test_per_variable_tolerance(self):
        write_sim(self.ref, {"elts/stress/stress.step1": REF_STRESS,
                             "elts/strain/strain.step1": REF_STRESS})
        write_sim(self.test, {"elts/stress/stress.step1": scaled(REF_STRESS, 1.001),
                              "elts/strain/strain.step1": scaled(REF_STRESS, 1.001)})
        opts = compare_sim.Options(per_var={"stress": compare_sim.Tolerance(1.0e-2, 1.0e-6)})
        n_fail, results = self.run_compare(opts)
        self.assertEqual(n_fail, 1)
        self.assertEqual([r.path for r in results if not r.passed], ["elts/strain/strain.step1"])

    def test_remove_mean_translation(self):
        """--remove-mean (step 8.3b): a uniform offset of a node variable is not a difference,
        but a non-uniform one still is; other variables are unaffected."""
        ref = [[0.0, 0.0, 1.0], [1.0, 0.5, 2.0], [2.0, 1.0, 0.5]]
        shifted = [[x + 0.3, y - 0.2, z + 0.1] for x, y, z in ref]
        write_sim(self.ref, {"nodes/coo/coo.step1": ref})
        write_sim(self.test, {"nodes/coo/coo.step1": shifted})
        self.assertEqual(self.run_compare()[0], 1)
        self.assertEqual(self.run_compare(compare_sim.Options(remove_mean=["coo"]))[0], 0)
        self.assertEqual(self.run_compare(compare_sim.Options(remove_mean=["vel"]))[0], 1)
        shifted[1][0] += 0.01
        write_sim(self.test, {"nodes/coo/coo.step1": shifted})
        self.assertEqual(self.run_compare(compare_sim.Options(remove_mean=["coo"]))[0], 1)

    def test_fortran_d_exponent_is_read(self):
        path = os.path.join(self.tmp.name, "d.txt")
        with open(path, "w", encoding="ascii") as handle:
            handle.write("0.1D+01 -.25d-01\n")
        self.assertEqual(compare_sim.read_values(path), [[1.0, -0.025]])

    def test_three_digit_exponent_is_read(self):
        path = os.path.join(self.tmp.name, "e3.txt")
        with open(path, "w", encoding="ascii") as handle:
            handle.write("-.6982454-309 0.1000000+101 1.5E-03\n")
        values = compare_sim.read_values(path)[0]
        self.assertAlmostEqual(values[0] / -0.6982454e-309, 1.0, places=12)
        self.assertEqual(values[1], 0.1e101)
        self.assertEqual(values[2], 1.5e-3)

    def test_comment_lines_compared_as_text(self):
        header = "% step     incr     Fx\n"
        for name, text in (("ref.sim", header), ("test.sim", header)):
            path = os.path.join(self.tmp.name, name, "results", "forces", "x0")
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w", encoding="ascii") as handle:
                handle.write(text + "  1  1  0.20084E+03\n")
        self.assertEqual(self.run_compare()[0], 0)
        path = os.path.join(self.test, "results", "forces", "x0")
        with open(path, "w", encoding="ascii") as handle:
            handle.write("% step     incr     Fy\n  1  1  0.20084E+03\n")
        n_fail, results = self.run_compare()
        self.assertEqual(n_fail, 1)
        self.assertIn("comment/text differs", results[0].message)

    def test_command_line_exit_status(self):
        write_sim(self.ref, {"elts/stress/stress.step1": REF_STRESS})
        write_sim(self.test, {"elts/stress/stress.step1": REF_STRESS})
        self.assertEqual(compare_sim.main([self.test, self.ref, "--quiet"]), 0)
        write_sim(self.test, {"elts/stress/stress.step1": scaled(REF_STRESS, 1.1)})
        self.assertEqual(compare_sim.main([self.test, self.ref, "--quiet"]), 1)
        self.assertEqual(compare_sim.main([self.test, os.path.join(self.tmp.name, "none")]), 2)
        self.assertEqual(compare_sim.main([self.test, self.ref, "--quiet",
                                           "--var-tol", "stress=0.2,1e-6"]), 0)


if __name__ == "__main__":
    unittest.main()
