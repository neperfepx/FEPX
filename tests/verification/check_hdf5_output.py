#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Check the HDF5 output of a run (output_format hdf5, step 9.1, ROADMAP D93).

Source : FEPX verification tests (ROADMAP step 9.1; layout: src/io/sim_writer_mod.f90 and
         dev/io_format.md §5). Development tool (CLAUDE.md §6): Python 3 + NumPy, the HDF5
         file read with the HDF5 command-line tools h5ls and h5dump (hdf5-tools).

Usage:
  check_hdf5_output.py SIMDIR --mesh MSH
      a run with output_format ascii hdf5: the HDF5 file against the text files of the same
      run, bit for bit, and against the mesh file:
        1. the root attributes: format, sizes, element type, the result lists and the steps
           of the text results;
        2. the mesh: /mesh/coo the node coordinates of MSH, /mesh/elts its 3D elements' nodes
           in its order, /mesh/elset their elset (first tag), /mesh/phase one phase per elset;
           /mesh/elts_xdmf the same nodes 0-based in the XDMF (VTK) order, each edge node at
           the middle of its VTK edge (straight-sided mesh, to 1e-12 of the size);
        3. every node and element variable at every step written: the dataset equals the
           text file (es24.16 is an exact round trip); a slip-system variable has the text
           line's values then NaN, and num_values the text line's counts;
        4. the tables of rows: each face set's forces, the convergence record (nr, attempt,
           accepted, bc_iter, solve columns mapped from the text columns), the macroscopic
           file; each with its "columns" attribute;
        5. /results/steps: the steps of the text files, at the time of the last force row of
           the step (when forces is printed; step 0 at time 0);
        6. the XDMF file (check_xdmf).
  check_hdf5_output.py SIMDIR --hdf5-only --reference REFDIR
      a run with output_format hdf5 alone: no .sim summary and no results/ directory; the HDF5
      file identical to REFDIR's (the same case with ascii hdf5; h5diff, every dataset and
      attribute) and the same XDMF file; the XDMF file checked (check_xdmf).

check_xdmf: an XDMF 3 file with a temporal collection, one uniform grid per row of
/results/steps with its time; topology Tet_10 (Hex_20) of num_elts elements on
/mesh/elts_xdmf; geometry the step's coo when it exists, else /mesh/coo; attributes elset,
phase and every printed variable written at the step; every data item an existing dataset of
the declared dimensions and number type.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ElementTree

import numpy as np

# XDMF (VTK) edge-node order: VTK edges (1-based corner nodes) of the quadratic tetrahedron
# and hexahedron, nodes 5-10 and 9-20 (VTK file formats; sim_writer_mod).
VTK_EDGES = {
    'tet10': [(1, 2), (2, 3), (3, 1), (1, 4), (2, 4), (3, 4)],
    'hex20': [(1, 2), (2, 3), (3, 4), (4, 1), (5, 6), (6, 7), (7, 8), (8, 5), (1, 5), (2, 6),
              (3, 7), (4, 8)],
}
NUM_CORNERS = {'tet10': 4, 'hex20': 8}
FORCE_COLUMNS = 'step incr Fx Fy Fz area time'
CONV_COLUMNS = ('incr iter nr r_norm rx_norm f_norm delu_norm delux_norm u_norm cg_iter '
                'attempt accepted bc_iter solve')
XDMF_VECTORS = {'rotrate', 'rotrate_spin', 'rotrate_slip'}

failures = []


def check(condition, message):
    """Record a failed check."""
    if not condition:
        failures.append(message)
    return condition


def run(command):
    """Run a command; its standard output (a failure ends the check)."""
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit('check_hdf5_output: %s failed:\n%s%s' % (' '.join(command), result.stdout,
                                                          result.stderr))
    return result.stdout


class H5File:
    """Read access to an HDF5 file through h5ls and h5dump."""

    def __init__(self, path):
        self.path = path
        self.datasets = {}
        # h5ls -r: "<path>  Dataset {rows, cols}" (an extendable table "{rows/Inf, cols}").
        for line in run(['h5ls', '-r', path]).splitlines():
            match = re.match(r'^(\S+)\s+Dataset \{(\d+)(?:/\w+)?, (\d+)(?:/\w+)?\}', line)
            if match:
                self.datasets[match.group(1)] = (int(match.group(2)), int(match.group(3)))
        self.tmp = tempfile.mkdtemp(prefix='check_hdf5_')

    def read(self, name):
        """A dataset as a (rows, cols) array (its file type: 64-bit float or 32-bit int)."""
        rows, cols = self.datasets[name]
        header = run(['h5dump', '-H', '-d', name, self.path])
        dtype = '<f8' if 'H5T_IEEE_F64LE' in header else '<i4'
        out = os.path.join(self.tmp, 'dataset.bin')
        if rows * cols == 0:
            return np.zeros((rows, cols), dtype=dtype)
        run(['h5dump', '-d', name, '-b', 'LE', '-o', out, self.path])
        values = np.fromfile(out, dtype=dtype)
        return values.reshape(rows, cols)

    def attribute(self, name):
        """An attribute's value: its text, or its integers."""
        text = run(['h5dump', '-a', name, self.path])
        data = text[text.index('DATA {'):]
        match = re.search(r'"(.*)"', data, re.S)
        if match:
            return match.group(1)
        return [int(v) for v in re.findall(r'(?<=[\s,(])-?\d+(?=[\s,])', data.split(':', 1)[-1])]


def read_text_table(path, skip=0):
    """The rows of a text file as lists of tokens (after skip lines)."""
    with open(path) as handle:
        return [line.split() for line in handle.readlines()[skip:] if line.strip()]


def same(a, b):
    """Two float arrays equal bit for bit (NaN equal to NaN)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return a.shape == b.shape and np.array_equal(a.view(np.int64), b.view(np.int64))


def read_mesh(path):
    """Node coordinates and the 3D elements (nodes, first tag) of a Neper msh file."""
    with open(path) as handle:
        lines = handle.read().splitlines()
    i = lines.index('$Nodes')
    num_nodes = int(lines[i + 1])
    coo = np.array([[float(v) for v in lines[i + 2 + j].split()[1:4]]
                    for j in range(num_nodes)])
    i = lines.index('$Elements')
    elts, elsets = [], []
    for j in range(int(lines[i + 1])):
        tokens = [int(v) for v in lines[i + 2 + j].split()]
        if tokens[1] in (11, 17):
            num_tags = tokens[2]
            elsets.append(tokens[3])
            elts.append(tokens[3 + num_tags:])
    return coo, np.array(elts), np.array(elsets)


def read_result_lists(sim_dir):
    """The results in the text output: the sorted lists of node and element results (forces and
    convergence are among the latter) and the largest step written."""
    results_dir = os.path.join(sim_dir, 'results')
    results = {'node': sorted(os.listdir(os.path.join(results_dir, 'nodes'))),
               'elt': sorted(os.listdir(os.path.join(results_dir, 'elts'))
                             + [name for name in ('forces', 'convergence')
                                if os.path.exists(os.path.join(results_dir, name))])}
    steps = [int(name.rsplit('.step', 1)[1])
             for entity in ('nodes', 'elts')
             for var in os.listdir(os.path.join(results_dir, entity))
             for name in os.listdir(os.path.join(results_dir, entity, var))]
    return results, max(steps)


def check_xdmf(sim_dir, h5):
    """The XDMF file against the HDF5 file."""
    path = os.path.join(sim_dir, 'simulation.xdmf')
    if not check(os.path.exists(path), 'no simulation.xdmf'):
        return
    root = ElementTree.parse(path).getroot()
    check(root.tag == 'Xdmf' and root.get('Version') == '3.0', 'XDMF root element')
    collection = root.find('Domain/Grid')
    check(collection.get('GridType') == 'Collection'
          and collection.get('CollectionType') == 'Temporal', 'temporal collection')
    steps = h5.read('/results/steps')
    grids = collection.findall('Grid')
    check(len(grids) == steps.shape[0], 'one XDMF grid per step written')
    element_type = h5.attribute('/element_type')
    num_elts = h5.attribute('/num_elts')[0]
    node_results = h5.attribute('/node_results').split()
    elt_results = [v for v in h5.attribute('/elt_results').split()
                   if v not in ('forces', 'convergence')]
    for grid, (step, time) in zip(grids, steps):
        step = int(step)
        label = 'XDMF step %d' % step
        check(grid.get('Name') == 'step%d' % step, label + ': grid name')
        check(float(grid.find('Time').get('Value')) == time, label + ': time')
        topology = grid.find('Topology')
        check(topology.get('TopologyType') == {'tet10': 'Tet_10', 'hex20': 'Hex_20'}[
            element_type], label + ': topology type')
        check(int(topology.get('NumberOfElements')) == num_elts, label + ': element count')
        items = {'topology': topology.find('DataItem'),
                 'geometry': grid.find('Geometry/DataItem')}
        check(items['topology'].text.endswith(':/mesh/elts_xdmf'), label + ': connectivity')
        coo = '/results/nodes/coo/coo.step%d' % step
        expect_geometry = coo if coo in h5.datasets else '/mesh/coo'
        check(items['geometry'].text == 'simulation.h5:' + expect_geometry, label + ': geometry')
        names = []
        for attribute in grid.findall('Attribute'):
            names.append(attribute.get('Name'))
            items['attribute ' + attribute.get('Name')] = attribute.find('DataItem')
            dims = h5.datasets.get(attribute.find('DataItem').text.split(':', 1)[1])
            name = attribute.get('Name')
            if name in node_results:
                kind, center = 'Vector', 'Node'
            elif dims[1] == 1:
                kind, center = 'Scalar', 'Cell'
            elif dims[1] == 3 and name in XDMF_VECTORS:
                kind, center = 'Vector', 'Cell'
            else:
                kind, center = 'Matrix', 'Cell'
            check(attribute.get('AttributeType') == kind and attribute.get('Center') == center,
                  label + ': attribute %s type' % name)
        expected = ['elset', 'phase']
        expected += [v for v in node_results
                     if '/results/nodes/%s/%s.step%d' % (v, v, step) in h5.datasets]
        expected += [v for v in elt_results
                     if '/results/elts/%s/%s.step%d' % (v, v, step) in h5.datasets]
        check(sorted(names) == sorted(expected), label + ': attributes %s, expected %s'
              % (names, expected))
        # Every data item: an existing dataset of the declared dimensions and type.
        for what, item in items.items():
            file_name, dataset = item.text.split(':', 1)
            check(file_name == 'simulation.h5' and item.get('Format') == 'HDF',
                  label + ': %s file' % what)
            if not check(dataset in h5.datasets, label + ': %s: no %s' % (what, dataset)):
                continue
            check(item.get('Dimensions') == '%d %d' % h5.datasets[dataset],
                  label + ': %s dimensions' % what)
            is_int = dataset.startswith('/mesh/') and dataset != '/mesh/coo'
            check(item.get('NumberType') == ('Int' if is_int else 'Float')
                  and item.get('Precision') == ('4' if is_int else '8'),
                  label + ': %s number type' % what)


def check_mesh(h5, mesh_path):
    """The mesh datasets against the mesh file."""
    coo, elts, elsets = read_mesh(mesh_path)
    check(same(h5.read('/mesh/coo'), coo), '/mesh/coo differs from the mesh file')
    stored = h5.read('/mesh/elts')
    check(np.array_equal(stored, elts), '/mesh/elts differs from the mesh file')
    check(np.array_equal(h5.read('/mesh/elset')[:, 0], elsets), '/mesh/elset')
    phase = h5.read('/mesh/phase')[:, 0]
    for elset in set(elsets):
        check(len(set(phase[elsets == elset])) == 1, '/mesh/phase: one phase per elset')
    # The XDMF connectivity: the corners, then each edge node at the middle of its VTK edge.
    element_type = h5.attribute('/element_type')
    xdmf = h5.read('/mesh/elts_xdmf') + 1
    corners = NUM_CORNERS[element_type]
    check(np.array_equal(xdmf[:, :corners], elts[:, :corners]), 'XDMF corner nodes')
    size = np.ptp(coo, axis=0).max()
    for k, (a, b) in enumerate(VTK_EDGES[element_type]):
        middle = 0.5 * (coo[xdmf[:, a - 1] - 1] + coo[xdmf[:, b - 1] - 1])
        error = np.abs(coo[xdmf[:, corners + k] - 1] - middle).max()
        check(error <= 1e-12 * size, 'XDMF edge node %d not at the middle of edge (%d, %d): %g'
              % (corners + k + 1, a, b, error))
    check(sorted(xdmf[0]) == sorted(elts[0]), 'XDMF nodes are a permutation of the element\'s')


def check_fields(sim_dir, h5, results):
    """Every node and element variable at every step against its text file."""
    steps_written = set()
    for entity, group in (('node', 'nodes'), ('elt', 'elts')):
        for var in results.get(entity, []):
            if var in ('forces', 'convergence'):
                continue
            directory = os.path.join(sim_dir, 'results', group, var)
            for file_name in sorted(os.listdir(directory)):
                step = int(file_name.rsplit('.step', 1)[1])
                steps_written.add(step)
                dataset = '/results/%s/%s/%s' % (group, var, file_name)
                if not check(dataset in h5.datasets, 'no dataset %s' % dataset):
                    continue
                rows = read_text_table(os.path.join(directory, file_name))
                table = h5.read(dataset)
                num_values = '/results/elts/%s/num_values' % var
                if entity == 'elt' and num_values in h5.datasets:
                    counts = h5.read(num_values)[:, 0]
                    check(list(counts) == [len(r) for r in rows],
                          '%s: num_values differs from the text lines' % dataset)
                    ok = all(same(table[i, :len(r)], [float(v) for v in r])
                             and np.isnan(table[i, len(r):]).all() for i, r in enumerate(rows))
                else:
                    ok = same(table, [[float(v) for v in r] for r in rows])
                check(ok, '%s differs from %s' % (dataset, file_name))
            # Nothing in the HDF5 file that the text files lack.
            prefix = '/results/%s/%s/' % (group, var)
            extra = [d for d in h5.datasets if d.startswith(prefix) and not d.endswith(
                'num_values') and not os.path.exists(os.path.join(directory,
                                                                  d[len(prefix):]))]
            check(not extra, 'datasets without a text file: %s' % extra)
    return steps_written


def check_tables(sim_dir, h5, results):
    """The tables of rows against the text files; the step times (from the forces)."""
    elt_results = results.get('elt', [])
    last_force_time = {0: 0.0}
    # Forces: v2's two header lines, then step incr Fx Fy Fz area time.
    if 'forces' in elt_results:
        directory = os.path.join(sim_dir, 'results', 'forces')
        for faset in sorted(os.listdir(directory)):
            dataset = '/results/forces/' + faset
            rows = [[float(v) for v in r] for r in
                    read_text_table(os.path.join(directory, faset), skip=2)]
            check(same(h5.read(dataset), rows), '%s differs from the text file' % dataset)
            check(h5.attribute(dataset + '/columns') == FORCE_COLUMNS, dataset + ' columns')
            for r in rows:
                last_force_time[int(r[0])] = r[6]
    # Convergence: the header line(s), then incr iter SA|NR six norms cg [attempt A|R
    # [bc_iter solve]].
    if 'convergence' in elt_results:
        expect = []
        for r in read_text_table(os.path.join(sim_dir, 'results', 'convergence')):
            if r[0].startswith('%'):
                continue
            row = [float(r[0]), float(r[1]), 1.0 if r[2] == 'NR' else 0.0]
            row += [float(v) for v in r[3:9]] + [float(r[9])]
            row += [float(r[10]), 1.0 if r[11] == 'A' else 0.0] if len(r) > 10 else [1.0, 1.0]
            row += [float(r[12]), float(r[13])] if len(r) > 12 else [0.0, 0.0]
            expect.append(row)
        check(same(h5.read('/results/convergence'), expect), '/results/convergence differs')
        check(h5.attribute('/results/convergence/columns') == CONV_COLUMNS,
              '/results/convergence columns')
    # Macroscopic file (periodic conditions): one header line.
    macro = os.path.join(sim_dir, 'results', 'macro')
    if os.path.exists(macro):
        rows = [[float(v) for v in r] for r in read_text_table(macro, skip=1)]
        check(same(h5.read('/results/macro'), rows), '/results/macro differs')
        check(len(h5.attribute('/results/macro/columns').split()) == 27, 'macro columns')
    else:
        check('/results/macro' not in h5.datasets, '/results/macro without results/macro')
    return last_force_time if 'forces' in elt_results else None


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('sim_dir')
    parser.add_argument('--mesh')
    parser.add_argument('--hdf5-only', action='store_true')
    parser.add_argument('--reference')
    args = parser.parse_args()
    h5 = H5File(os.path.join(args.sim_dir, 'simulation.h5'))

    if args.hdf5_only:
        # No text output; the same HDF5 and XDMF files as the reference run's.
        check(not os.path.exists(os.path.join(args.sim_dir, 'results')), 'results/ exists')
        check(os.path.isdir(os.path.join(args.sim_dir, 'inputs')), 'no inputs/')
        result = subprocess.run(['h5diff', os.path.join(args.reference, 'simulation.h5'),
                                 h5.path], capture_output=True, text=True)
        check(result.returncode == 0, 'the HDF5 file differs from the reference run\'s:\n'
              + result.stdout + result.stderr)
        with open(os.path.join(args.reference, 'simulation.xdmf')) as a, \
                open(os.path.join(args.sim_dir, 'simulation.xdmf')) as b:
            check(a.read() == b.read(), 'the XDMF file differs from the reference run\'s')
    else:
        # 1. The summary.
        results, num_steps = read_result_lists(args.sim_dir)
        check(h5.attribute('/format') == 'fepx-hdf5 1', 'format attribute')
        check(sorted(h5.attribute('/node_results').split()) == results['node'], 'node_results')
        check(sorted(h5.attribute('/elt_results').split()) == results['elt'], 'elt_results')
        check(h5.attribute('/num_steps')[0] == num_steps, 'num_steps')
        check(h5.attribute('/num_nodes')[0] == h5.datasets['/mesh/coo'][0], 'num_nodes')
        check(h5.attribute('/num_elts')[0] == h5.datasets['/mesh/elts'][0], 'num_elts')
        # 2.-4. The mesh, the fields, the tables of rows.
        check_mesh(h5, args.mesh)
        steps_written = check_fields(args.sim_dir, h5, results)
        step_times = check_tables(args.sim_dir, h5, results)
        # 5. The steps written and their times.
        steps = h5.read('/results/steps')
        check(sorted(steps_written) == list(steps[:, 0].astype(int)), '/results/steps: %s, '
              'text files %s' % (list(steps[:, 0]), sorted(steps_written)))
        if step_times is not None:
            for step, time in steps:
                check(time == step_times[int(step)], '/results/steps: time of step %d' % step)
        check(h5.attribute('/results/steps/columns') == 'step time', '/results/steps columns')
    # 6. The XDMF file.
    check_xdmf(args.sim_dir, h5)

    for failure in failures:
        print('FAIL ' + failure)
    print('check_hdf5_output: %d datasets, %s' % (len(h5.datasets),
                                                  'FAILED' if failures else 'passed'))
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
