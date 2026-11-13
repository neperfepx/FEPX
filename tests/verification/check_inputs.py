#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check the copies of the input files in simulation.sim/inputs.

Source : FEPX verification tests (the input copies of src/io/sim_writer_mod.f90, sim_open; the
input list of src/app/fepx.f90). Development tool (CLAUDE.md §6).

The directory SIM_DIR/inputs must hold exactly the files given by --expect, each identical to the
file of the same name in the run directory (the parent of SIM_DIR).

Usage: check_inputs.py SIM_DIR --expect FILE1 FILE2 ...
Exit status 0 if the check passes, 1 otherwise.
"""

import argparse
import filecmp
import os
import sys


def main():
    # Arguments: the simulation directory and the expected input files.
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('sim_dir')
    parser.add_argument('--expect', nargs='+', required=True)
    args = parser.parse_args()
    inputs_dir = os.path.join(args.sim_dir, 'inputs')
    run_dir = os.path.dirname(os.path.abspath(args.sim_dir))
    errors = []

    # The set of copied files.
    found = sorted(os.listdir(inputs_dir)) if os.path.isdir(inputs_dir) else []
    if found != sorted(args.expect):
        errors.append(f'inputs/ holds {found}, expected {sorted(args.expect)}')

    # Each copy identical to the run directory's file.
    for name in args.expect:
        copy = os.path.join(inputs_dir, name)
        if os.path.isfile(copy) and not filecmp.cmp(copy, os.path.join(run_dir, name),
                                                    shallow=False):
            errors.append(f'inputs/{name} differs from the run directory\'s {name}')

    for error in errors:
        print('FAIL:', error)
    if not errors:
        print('PASS: inputs/ holds', ', '.join(sorted(args.expect)))
    return 1 if errors else 0


if __name__ == '__main__':
    sys.exit(main())
