#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Helper of tutorials/tools/make_tutorial_pages.sh: the text and the curves of the manual's tutorial pages.

Step 9.11 (ROADMAP D95). A development tool (CLAUDE.md §6: Python 3 + NumPy, here with
Matplotlib for the curves); it computes nothing new: every curve is a column of a file FEPX
wrote (results/forces/<faset>, results/macro, results/elts/<var>) or of the volume averages that
Neper writes (`neper -S -resmesh stress,strain`: results/mesh/<var>/<var>.step<n>).

Subcommands:
  excerpt <log> <out>          the terminal output of a run, shortened as running.rst shows it:
                               from the banner to the end of the first increment, "[...]", then
                               from "Final step terminated" to the end
  tree <sim_dir> <out>         the tree of a simulation directory (the file format of
                               fileformat.rst), each result variable's step files shortened
  scale <result_file> <column> a color scale "<low>:<high>" for neper -V -dataeltscale
  plot <case> <out_dir> <sim_dir>...
                               the curves of a tutorial page, as PNG files in <out_dir>
  postprocess-test <reference_dir> <mesh> <script> <work_dir>
                               the docs test of the post-processing tutorial's script
  check <tutorials_txt> <pages_dir>
                               the consistency check of the docs test: every generated file a
                               page names exists, and every case of the list has its files
"""

import os
import re
import sys

import numpy as np

# Width of a terminal excerpt line (FEPX writes at most 72 columns, D94).
EXCERPT_GAP = "  [...]"

# Colors of the curves: the first slots of a validated categorical palette (blue, orange,
# aqua), assigned in this fixed order; each series also has its own marker.
SERIES_COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]
SERIES_MARKERS = ["o", "s", "^"]

# Stress and strain components of FEPX's 6-vectors (order 11 22 33 23 13 12, conventions.md).
COMPONENT_INDEX = {"11": 0, "22": 1, "33": 2, "23": 3, "13": 4, "12": 5}


def write_excerpt(log_path, out_path):
    """Write the shortened terminal output of a run.

    log_path  the run's terminal output (str, path)
    out_path  the excerpt written (str, path)
    Local variables: lines (list of str, the log's lines), increment_lines, step_lines (lists of
    int, indices of the lines that start an increment, a step), first_cut (int, index of the line
    that starts the second increment or the second step, whichever comes first), final_start
    (int, index of the "Final step terminated" line).
    """
    # Read the log, without trailing spaces, keeping FEPX's lines only: every line FEPX
    # writes starts with "Info   :", "Warning:", "Error  :" or is a rule of "=" (D94); a log
    # that also caught the MPI library's own messages (standard error) is thus cleaned.
    fepx_line = re.compile(r"^(Info   :|Warning:|Error  :|={20,})")
    with open(log_path) as log_file:
        lines = [line.rstrip() for line in log_file if fepx_line.match(line)]
    # The first increment ends where the run's second increment or second step starts (a
    # restarted run numbers them on from the first run).
    increment_lines = [index for index, line in enumerate(lines)
                       if re.match(r"^Info   :   - Increment \d+:", line)]
    step_lines = [index for index, line in enumerate(lines)
                  if re.match(r"^Info   : Running step \d+\.\.\.", line)]
    cuts = increment_lines[1:2] + step_lines[1:2]
    first_cut = min(cuts) if cuts else None
    # The closing summary starts at "Final step terminated".
    final_start = None
    for index, line in enumerate(lines):
        if "Final step terminated" in line:
            final_start = index
            break
    # Keep the whole log when it has no second increment (nothing to shorten).
    if first_cut is None or final_start is None or final_start <= first_cut:
        kept = lines
    else:
        kept = lines[:first_cut] + [EXCERPT_GAP] + lines[final_start:]
    with open(out_path, "w") as out_file:
        out_file.write("\n".join(kept) + "\n")


def tree_lines(path, prefix):
    """Return the lines of the tree under a directory (recursive).

    path    the directory (str)
    prefix  the text before each entry of this level (str)
    The step files of a result variable (<var>.step<n>) are shortened to the first two and the
    last when there are more than four. Hidden files (the .sim summary) are listed last.
    Local variables: names (list of str, the entries, sorted), step_names (list of str, the
    entries named <var>.step<n>, in step order), shown (list of str or None, the entries
    written, None for the "..." line).
    """
    names = sorted(os.listdir(path), key=lambda name: (name.startswith("."), name))
    # Order the step files by step number and shorten long lists.
    step_pattern = re.compile(r"^(.*)\.step(\d+)$")
    step_names = [name for name in names if step_pattern.match(name)]
    if step_names and len(step_names) == len(names):
        step_names.sort(key=lambda name: int(step_pattern.match(name).group(2)))
        shown = step_names if len(step_names) <= 4 else \
            step_names[:2] + [None] + step_names[-1:]
    else:
        shown = names
    # Write each entry, and the entries of each subdirectory under it.
    lines = []
    for index, name in enumerate(shown):
        last = index == len(shown) - 1
        branch = "`-- " if last else "|-- "
        lines.append(prefix + branch + ("..." if name is None else name))
        if name is not None and os.path.isdir(os.path.join(path, name)):
            lines += tree_lines(os.path.join(path, name), prefix + ("    " if last else "|   "))
    return lines


def write_tree(sim_dir, out_path):
    """Write the tree of a simulation directory (named simulation.sim on the page).

    sim_dir   the simulation directory (str)
    out_path  the tree written (str, path)
    """
    lines = ["simulation.sim"] + tree_lines(sim_dir, "")
    with open(out_path, "w") as out_file:
        out_file.write("\n".join(lines) + "\n")


def read_forces(sim_dir, faset):
    """Return the force table of a face set: columns step incr Fx Fy Fz area time.

    sim_dir  the simulation directory (str)
    faset    the face set name, e.g. "z1" (str)
    Returns a 2D float array, one row per increment (row 0: the initial state).
    """
    return np.atleast_2d(np.loadtxt(os.path.join(sim_dir, "results", "forces", faset),
                                    comments="%"))


def read_mesh_average(sim_dir, var):
    """Return Neper's volume averages of a 6-vector result at the printed steps.

    sim_dir  the simulation directory, post-processed by neper -S -resmesh (str)
    var      "stress" or "strain" (str)
    Returns (steps, values): steps (int array), values (float array, one row of 6 components,
    order 11 22 33 23 13 12, per step).
    """
    directory = os.path.join(sim_dir, "results", "mesh", var)
    steps = sorted(int(name.rsplit("step", 1)[1]) for name in os.listdir(directory))
    values = np.array([np.loadtxt(os.path.join(directory, f"{var}.step{step}"))
                       for step in steps])
    return np.array(steps), values


def face_stress(forces, axis):
    """Return the normal stress of a face, force over current area [MPa], per increment.

    forces  a force table of read_forces (2D float array)
    axis    0, 1 or 2 for the x, y or z force component (int)
    The initial row has a zero force; its stress is 0.
    """
    return forces[:, 2 + axis] / forces[:, 5]


def new_figure(num_panels):
    """Return a Matplotlib figure and its axes, in the pages' style.

    num_panels  1 or 2 panels side by side (int)
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.color": "#dddddd", "grid.linewidth": 0.6,
                         "lines.linewidth": 2.0, "lines.markersize": 5,
                         "legend.frameon": False})
    figure, axes = plt.subplots(1, num_panels, figsize=(4.2 * num_panels, 3.4),
                                constrained_layout=True)
    return figure, np.atleast_1d(axes)


def save_figure(figure, out_dir, name):
    """Write a figure as <out_dir>/<name>.png (150 dpi) and close it."""
    import matplotlib.pyplot as plt
    figure.savefig(os.path.join(out_dir, name + ".png"), dpi=150)
    plt.close(figure)


def plot_series(axis, x_values, y_values, slot, label=None, step_ends=False):
    """Draw one series in the color and marker of its slot.

    axis      the Matplotlib axes
    x_values, y_values  the series (1D float arrays)
    slot      the series' slot, 0, 1 or 2 (fixed order: color and marker) (int)
    label     the legend label (str or None)
    step_ends values at the ends of the steps only (bool): markers joined by a thin dotted line,
              which is not the path between them; otherwise a solid line through every value
    """
    if step_ends:
        axis.plot(x_values, y_values, color=SERIES_COLORS[slot], marker=SERIES_MARKERS[slot],
                  markersize=7, linestyle=":", linewidth=1.2, label=label)
    else:
        axis.plot(x_values, y_values, color=SERIES_COLORS[slot], marker=SERIES_MARKERS[slot],
                  label=label)


def plot_resmesh_stress_strain(axis, sim_dir, components, labels=None):
    """Draw the volume-averaged stress against strain of the given components.

    axis        the Matplotlib axes
    sim_dir     the simulation directory, post-processed by neper -S -resmesh (str)
    components  the components, e.g. ["33"] (list of str)
    labels      the legend labels (list of str, or None for the components)
    """
    _, stress = read_mesh_average(sim_dir, "stress")
    _, strain = read_mesh_average(sim_dir, "strain")
    for slot, component in enumerate(components):
        index = COMPONENT_INDEX[component]
        label = labels[slot] if labels else component
        plot_series(axis, 100.0 * strain[:, index], stress[:, index], slot, label,
                    step_ends=True)
    # One component: named on the axes (no legend); several: by the legend.
    name = f" {components[0]}" if len(components) == 1 else ""
    axis.set_xlabel(f"Strain{name} (volume average) [%]")
    axis.set_ylabel(f"Stress{name} (volume average) [MPa]")


def plot_face_stress_time(axis, sim_dir, fasets, labels=None):
    """Draw the normal stress of faces (force over area) against time, per increment.

    axis     the Matplotlib axes
    sim_dir  the simulation directory (str)
    fasets   the face sets, e.g. ["x1", "y1", "z1"] (list of str); the normal of xN is x, etc.
    labels   the legend labels (list of str, or None for the face sets)
    """
    for slot, faset in enumerate(fasets):
        forces = read_forces(sim_dir, faset)
        axis_index = "xyz".index(faset[0])
        label = labels[slot] if labels else faset
        plot_series(axis, forces[:, 6], face_stress(forces, axis_index), slot, label)
    axis.set_xlabel("Time [s]")
    if len(fasets) == 1:
        axis.set_ylabel(f"Force {fasets[0][0]} / area of {fasets[0]} [MPa]")
    else:
        axis.set_ylabel("Normal force / area of the face [MPa]")


def plot_case(case, out_dir, sim_dirs):
    """Write the curves of a tutorial page.

    case      the case name of tutorials.txt, or a page's curve name (str)
    out_dir   the directory of the PNG files (str)
    sim_dirs  the simulation directories the curves read (list of str), per case:
              simple_rve: simple_simulation, simple_rve, simple_rve_grip; hcp: basal,
              prismatic, pyramidal; others: the case's own
    """
    os.makedirs(out_dir, exist_ok=True)
    sim_dir = sim_dirs[0]
    if case in ("simple_simulation", "mixed_targets", "cruciform"):
        # Stress-strain from the volume averages, and the face stress per increment.
        figure, axes = new_figure(2)
        if case == "cruciform":
            plot_resmesh_stress_strain(axes[0], sim_dir, ["11", "22"])
            plot_face_stress_time(axes[1], sim_dir, ["x1", "y1"])
            for axis in axes:
                axis.legend()
        else:
            plot_resmesh_stress_strain(axes[0], sim_dir, ["33"])
            plot_face_stress_time(axes[1], sim_dir, ["z1"])
        save_figure(figure, out_dir, "curves")
    elif case == "simple_rve":
        # The three uniaxial inputs on the same axes.
        labels = ["general (simple simulation)", "uniaxial_minimal", "uniaxial_grip"]
        figure, axes = new_figure(2)
        for slot, (directory, label) in enumerate(zip(sim_dirs, labels)):
            _, stress = read_mesh_average(directory, "stress")
            _, strain = read_mesh_average(directory, "strain")
            plot_series(axes[0], 100.0 * strain[:, 2], stress[:, 2], slot, label,
                        step_ends=True)
            forces = read_forces(directory, "z1")
            plot_series(axes[1], forces[:, 6], face_stress(forces, 2), slot, label)
        axes[0].set_xlabel("Strain 33 (volume average) [%]")
        axes[0].set_ylabel("Stress 33 (volume average) [MPa]")
        axes[1].set_xlabel("Time [s]")
        axes[1].set_ylabel("Force z / area of z1 [MPa]")
        axes[0].legend()
        save_figure(figure, out_dir, "curves")
    elif case == "periodic_rve":
        # The macroscopic stress 33 against the macroscopic strain Fbar_33 - 1, per increment.
        macro = np.atleast_2d(np.loadtxt(os.path.join(sim_dir, "results", "macro"),
                                         comments="%"))
        figure, axes = new_figure(1)
        strain_33 = macro[:, 3 + 8] - 1.0
        stress = macro[:, 3 + 9 + 9:]
        for slot, component in enumerate(["11", "22", "33"]):
            plot_series(axes[0], 100.0 * strain_33, stress[:, COMPONENT_INDEX[component]],
                        slot, component)
        axes[0].set_xlabel("Macroscopic strain 33 [%]")
        axes[0].set_ylabel("Macroscopic stress [MPa]")
        axes[0].legend()
        save_figure(figure, out_dir, "curves")
    elif case in ("multiaxial_csr", "multiaxial_clr", "restart"):
        # The three normal components: stress-strain at the steps, face stress per increment.
        figure, axes = new_figure(2)
        plot_resmesh_stress_strain(axes[0], sim_dir, ["11", "22", "33"])
        plot_face_stress_time(axes[1], sim_dir, ["x1", "y1", "z1"])
        for axis in axes:
            axis.legend()
        save_figure(figure, out_dir, "curves")
    elif case == "overwriting":
        figure, axes = new_figure(1)
        plot_face_stress_time(axes[0], sim_dir, ["z1"])
        save_figure(figure, out_dir, "curves")
    elif case == "hcp":
        # Stress-strain of the three orientations, and the slip rate of each family.
        labels = ["basal", "prismatic", "pyramidal <c+a>"]
        families = [range(0, 3), range(3, 6), range(6, 18)]
        figure, axes = new_figure(2)
        plot_resmesh_stress_strain(axes[0], sim_dirs[0], ["33"], [labels[0]])
        for slot in (1, 2):
            _, stress = read_mesh_average(sim_dirs[slot], "stress")
            _, strain = read_mesh_average(sim_dirs[slot], "strain")
            plot_series(axes[0], 100.0 * strain[:, 2], stress[:, 2], slot, labels[slot],
                        step_ends=True)
        axes[0].legend(loc="lower right")
        # Mean over the elements of each family's summed slip-rate magnitudes, per orientation:
        # the family the orientation favors only (the others stay at 0).
        for slot, directory in enumerate(sim_dirs):
            steps, strain = read_mesh_average(directory, "strain")
            rates = []
            for step in steps:
                sliprate = np.abs(np.loadtxt(os.path.join(
                    directory, "results", "elts", "sliprate", f"sliprate.step{step}")))
                rates.append(sliprate[:, list(families[slot])].sum(axis=1).mean())
            plot_series(axes[1], 100.0 * strain[:, 2], rates, slot, labels[slot],
                        step_ends=True)
        axes[1].set_xlabel("Strain 33 (volume average) [%]")
        axes[1].set_ylabel("Slip rate of the favored family [1/s]")
        axes[1].legend(loc="lower right")
        save_figure(figure, out_dir, "curves")
    else:
        raise SystemExit(f"tutorial_pages.py: no curves for case {case}")


def color_scale(result_file, column):
    """Return a color scale "<t0>:<t1>:...:<tn>" for a field image (neper -V -dataeltscale).

    result_file  an element result file, results/elts/<var>/<var>.step<n> (str)
    column       the 1-based column (int)
    The scale spans the 2nd to 98th percentiles of the column (so that a few extreme elements
    do not wash out the field), widened to round ticks: a step of 1, 2, 2.5 or 5 times a power of
    ten giving 4 to 8 intervals; a scale that crosses zero is made symmetric. The ticks are
    written with the decimals of the step, which Neper keeps on the scale's labels.
    Local variables: values (1D float array), low, high (float, the percentiles and then the end
    ticks), step (float, the tick step), decimals (int, decimals of the labels).
    """
    values = np.loadtxt(result_file, ndmin=2)[:, column - 1]
    low, high = np.percentile(values, [2.0, 98.0])
    if low < 0.0 < high:
        high = max(-low, high)
        low = -high
    if high <= low:
        high = low + (abs(low) if low != 0.0 else 1.0)
    # The smallest round step that gives at most 8 intervals.
    magnitude = 10.0 ** np.floor(np.log10((high - low) / 8.0))
    for factor in (1.0, 2.0, 2.5, 5.0, 10.0):
        step = factor * magnitude
        if (np.ceil(high / step - 1e-9) - np.floor(low / step + 1e-9)) <= 8:
            break
    low = np.floor(low / step + 1e-9) * step
    high = np.ceil(high / step - 1e-9) * step
    num_intervals = int(round((high - low) / step))
    # Decimals of the step (2.5 x 10^k needs one more than 10^k).
    decimals = max(0, int(-np.floor(np.log10(step) + 1e-9)) + (1 if factor == 2.5 else 0))
    ticks = [low + index * step for index in range(num_intervals + 1)]
    return ":".join(f"{tick:.{decimals}f}" for tick in ticks)


def postprocess_test(reference_dir, mesh_file, script, work_dir):
    """Run the post-processing tutorial's script on reference results (the docs test).

    reference_dir  a reference of tests/regression/reference/fepx (its results/) (str)
    mesh_file      the case's mesh, copied as simulation.sim/inputs/simulation.msh (str)
    script         tutorials/postprocessing/postprocess.py (str)
    work_dir       a scratch directory, replaced (str)
    Returns 0 when the script ran and wrote its two figures, 1 otherwise.
    """
    import shutil
    import subprocess
    # A simulation directory made of the reference results and the mesh.
    shutil.rmtree(work_dir, ignore_errors=True)
    sim_dir = os.path.join(work_dir, "simulation.sim")
    shutil.copytree(os.path.join(reference_dir, "results"), os.path.join(sim_dir, "results"))
    os.makedirs(os.path.join(sim_dir, "inputs"))
    shutil.copy(mesh_file, os.path.join(sim_dir, "inputs", "simulation.msh"))
    # Run the script there, as the page says.
    status = subprocess.run([sys.executable, os.path.abspath(script), "simulation.sim"],
                            cwd=work_dir).returncode
    missing = [name for name in ("stress_time.png", "grains.png")
               if not os.path.isfile(os.path.join(work_dir, name))]
    if status != 0 or missing:
        print(f"postprocess.py: status {status}, missing {missing}")
        return 1
    return 0


def check(tutorials_txt, pages_dir):
    """Check that the pages' generated files exist (the docs test; no Neper needed).

    tutorials_txt  the case list (str, path)
    pages_dir      doc (str, path)
    Every "[...]generated/..." path named in an .rst file of the manual must exist (relative to
    the page), and every case of the list must have generated/<case>/terminal.txt and tree.txt.
    Returns the number of problems (int); each is printed.
    """
    problems = 0
    generated_dir = os.path.join(pages_dir, "tutorials", "generated")
    # Every case of the list has its terminal output and tree.
    with open(tutorials_txt) as list_file:
        cases = [line.split()[0] for line in list_file
                 if line.strip() and not line.startswith("#")]
    for case in cases:
        for name in ("terminal.txt", "tree.txt"):
            if not os.path.isfile(os.path.join(generated_dir, case, name)):
                print(f"missing: tutorials/generated/{case}/{name}")
                problems += 1
    # Every generated file a page names exists.
    # (The path as written, relative to the page: "generated/..." in a tutorial page,
    # "tutorials/generated/..." in a page of the manual's top directory.)
    reference_pattern = re.compile(r"([A-Za-z0-9_./-]*generated/[A-Za-z0-9_./-]+\.(?:png|txt))")
    for root, _, files in os.walk(pages_dir):
        for name in files:
            if not name.endswith(".rst"):
                continue
            path = os.path.join(root, name)
            with open(path) as rst_file:
                for target in reference_pattern.findall(rst_file.read()):
                    if not os.path.isfile(os.path.join(root, target)):
                        print(f"{os.path.relpath(path, pages_dir)}: missing {target}")
                        problems += 1
    return problems


def main(arguments):
    """Run a subcommand (see the module description)."""
    if len(arguments) >= 3 and arguments[0] == "excerpt":
        write_excerpt(arguments[1], arguments[2])
    elif len(arguments) >= 3 and arguments[0] == "tree":
        write_tree(arguments[1], arguments[2])
    elif len(arguments) >= 4 and arguments[0] == "plot":
        plot_case(arguments[1], arguments[2], arguments[3:])
    elif len(arguments) == 3 and arguments[0] == "scale":
        print(color_scale(arguments[1], int(arguments[2])))
    elif len(arguments) == 5 and arguments[0] == "postprocess-test":
        return postprocess_test(*arguments[1:])
    elif len(arguments) == 3 and arguments[0] == "check":
        problems = check(arguments[1], arguments[2])
        print(f"tutorial pages: {problems} problem(s)")
        return 1 if problems else 0
    else:
        print(__doc__)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
