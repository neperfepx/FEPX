.. _special_commands:

Special Commands (:file:`simulation.cfg`)
=========================================

Two kinds of commands change where a simulation takes its input from: the input files that
replace the orientations, phases or parameters of the mesh file and the configuration, and the
restart of a simulation from the restart point of an earlier run.

.. _input_files:

External Input Files
--------------------

The orientations and phases are normally those of the mesh file. They can be replaced by those
of separate files (see :ref:`file_formats`):

.. confval:: read_ori_from_file

   The orientations are read from :file:`simulation.ori` instead of the mesh file.

   .. include:: generated/keywords/read_ori_from_file.inc

.. confval:: read_phase_from_file

   The phases are read from :file:`simulation.phase` instead of the mesh file.

   .. include:: generated/keywords/read_phase_from_file.inc

.. confval:: read_from_opt_file

   Per-element-set or per-element strengths and hardening parameters are read from
   :file:`simulation.opt` (:ref:`per_element_values`; format: :ref:`opt_file`).

   .. include:: generated/keywords/read_from_opt_file.inc

.. _restart:

Simulation Restart
------------------

With ``print restart``, FEPX writes a *restart point* at the end of every step (printed or
not): one file per process, :file:`simulation.sim/restart/restart.rank<r>`, holding everything
needed to continue (coordinates, velocities, the whole material state of every point, the
counters, the step just completed, the values the load-control drivers carry from step to
step, and a fingerprint of the configuration). Only the latest restart point is kept; it is
replaced safely (a run interrupted while writing leaves the previous one). A run that stops
within a step (a run limit, necking, a failed increment, a crash) leaves the restart point of
its last completed step.

To continue, write a new :file:`simulation.cfg` that lists **only the steps that follow** the
restart point, add ``restart on``, and run FEPX again in the same directory, on the same number
of processes:

.. confval:: restart

   ``on``: the run continues the simulation in :file:`simulation.sim` from its restart point
   (written by ``print restart``); the configuration lists only the steps that follow it.

   .. include:: generated/keywords/restart.inc


.. code-block:: none

  restart on
  ...                        (material, boundary conditions: unchanged)
  number_of_steps 2          (the steps after the restart point)
  strain_rate   33 -0.01 0.01
  target_strain 33 -0.004 0.004
  ...
  print restart              (to leave a restart point again)

Its step 1 continues the restart point's step N and is written as step N + 1; targets are
absolute (cumulative strain, load, stress, time), as in any configuration, and each step is
planned from where the previous one ended, so the loading may continue, reverse or switch to
another control mode. A test cycled in stages thus needs only the steps of the next stage.
The output continues in :file:`simulation.sim`: the step files are numbered on, and the force,
convergence and macroscopic files continue from their rows at the restart point (rows written
by an interrupted run after it are removed); with the HDF5 output the same holds in
:file:`simulation.h5` (see :ref:`output_formats`).

Rules:

- The mesh, the material, the boundary conditions (``boundary_conditions``, ``set_bc``), the
  constraints, the control component, the periodic kinematics, the input-file switches and the
  output format (``output_format``) must be those of the run that wrote the restart point: a
  difference is an error naming the line. For the output format only the formats written
  count: leaving ``output_format`` out and giving ``output_format ascii`` are the same.
  Other changes (solver and point-solve options, sub-incrementation, run limits, output
  requests) are allowed with a warning.
- The same number of processes (an error otherwise).
- ``restart on`` without a restart point is an error; without ``restart on`` an existing
  :file:`simulation.sim` is an error, as for any run.
- ``restart on`` does not imply ``print restart``: give both to leave a new restart point.

On the same number of processes, a run interrupted and continued this way gives the same
results, bit for bit, as the uninterrupted run, with two exceptions that depend on the whole
configuration: the columns of :file:`results/convergence` (a later part that needs more columns
writes a new header line), and, at a load or stress rate, the perturbation floor of the
boundary-velocity iteration (the largest target of each component; the restart point carries
that of the steps before it), which agree when the largest targets come first.
