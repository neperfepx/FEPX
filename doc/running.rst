.. _running:

Running a Simulation
====================

From the directory that holds :file:`simulation.cfg` and :file:`simulation.msh` (and the
:file:`simulation.ori`, :file:`simulation.phase` or :file:`simulation.opt` files that the
configuration asks for), run the command below. A tessellation :file:`simulation.tess` in the
directory is not read, but it is copied with the input files into :file:`simulation.sim/inputs`
(:ref:`simulation_directory`), for Neper's post-processing and visualization.

.. code-block:: console

  $ mpirun -np <N> fepx

FEPX must always be started through the MPI launcher (``mpirun`` or ``mpiexec``), also on one
process. On a cluster, use the launcher of the job scheduler. The only exceptions are the two
options below, which do not start MPI and are given without the launcher:

.. code-block:: console

  $ fepx --version
  3.0.1-dev (git 465464d2, 2026-10-26)
  $ fepx --help

``--version`` writes the version. A release writes it alone, e.g. ``3.0.0``. Otherwise, it is
followed by the git commit the program was built from (its short hash, with ``-dirty`` when the
source tree had uncommitted changes, and its date as YYYY-MM-DD): the development version, e.g.
``3.0.1-dev (git 465464d2, 2026-10-26)``, or, between releases on the main branch, the release
with the number of commits since it, e.g. ``3.0.0+3 (git 9f8e7d6c, 2026-11-20)``. A build
outside a git repository writes the version alone. ``--help`` writes the banner and the usage.
Any other argument is an error.

Number of Processes and Partitions
----------------------------------

The mesh is distributed as Neper partitioned it: process :math:`i` takes partition :math:`i`, so
the mesh must have as many partitions as there are processes. Partition a mesh for ``<N>``
processes when it is generated (``neper -M simulation.tess ... -part <N>``), or partition an
existing mesh with:

.. code-block:: console

  $ neper -M -loadmesh simulation.msh -part <N> -o simulation

(For a periodic mesh, this needs Neper 5.0.1 or later.)

If the number of processes differs from the number of partitions, or the mesh is not partitioned,
FEPX stops with an error that gives this command. The exception is a run on 1 process, which takes
any mesh whole; for a partitioned mesh FEPX warns:

.. code-block:: none

  Warning:     > The mesh has 8 partitions, but 1 rank was started: the
  Warning:       rank takes the whole mesh.

The results do not depend on the number of processes (to the solver tolerance).

.. _running_cluster:

Running on a Cluster
--------------------

On a cluster, FEPX is submitted to the job scheduler as a batch job. The scheduler's options
are too many to describe here; the script below, :file:`fepx_slurm.sh`, is an example for the
Slurm scheduler that runs FEPX on 8 nodes with 64 processes per node (512 processes), with
Open MPI (:file:`fepx` on the ``PATH`` of the compute nodes; otherwise give its full path). The
measurements of :ref:`performance_cluster` used a script of this form, with the Open MPI settings
given after the options:

.. code-block:: bash

  #!/bin/bash
  #SBATCH -J fepxjob
  #SBATCH -e error.%A
  #SBATCH -o output.%A
  #SBATCH --nodes=8
  #SBATCH --ntasks-per-node=64
  #SBATCH --exclusive

  echo "$SLURM_JOB_NODELIST are the nodes assigned"
  echo "$SLURM_NTASKS processes on $SLURM_NNODES nodes with $SLURM_NTASKS_PER_NODE per node"

  mpirun --bind-to core -np $SLURM_NTASKS fepx

  exit 0

It is submitted from the simulation directory with:

.. code-block:: console

  $ sbatch fepx_slurm.sh

The terminal output (:ref:`terminal_output`) goes to :file:`output.<job id>` and the error
messages to :file:`error.<job id>`. Partition the mesh for the number of processes of the job
(here 512), as described above.

What the options do, and what is known about them:

- ``--nodes`` and ``--ntasks-per-node``: use one process per physical core, not per hyperthread.
- ``--exclusive``: no other job on the nodes, which share memory bandwidth and the network card
  with whatever else runs there.
- ``--bind-to core``: holds each process on one core. It is good practice, but not a measured
  gain: on the cluster below, the processes were pinned to one core each without it as well.
- ``--constraint`` (add it if your cluster defines node features): keep a job on one type of
  processor. With Open MPI 5, :command:`mpirun` itself crashed at startup on some jobs that mixed
  two processor types.

**The settings of Open MPI.** ``--bind-to core`` is the only option on the command line above, but
the measurements also used settings from the system-wide file
:file:`/etc/openmpi/openmpi-mca-params.conf`, identical on every node (Open MPI 5.0, 10 Gb Ethernet
and no InfiniBand):

.. code-block:: text

  pml = ob1
  btl = self,sm,tcp
  btl_tcp_if_include = eno1np0
  btl_base_warn_component_unused = 0

.. list-table::
  :header-rows: 1
  :widths: 34 40 26

  * - Setting
    - What it does
    - Applies elsewhere?
  * - ``pml = ob1``
    - Selects Open MPI's point-to-point layer that runs over the transports named by ``btl``.
    - On clusters with Ethernet only. Not on InfiniBand or other fast fabrics, where it can keep
      Open MPI from using their usual layers (such as UCX).
  * - ``btl = self,sm,tcp``
    - The transports: a process to itself (``self``), shared memory within a node (``sm``) and TCP
      between nodes (``tcp``).
    - On clusters with Ethernet only; not on fast fabrics, since it leaves them unused. Open MPI 4
      calls shared memory ``vader``, not ``sm``.
  * - ``btl_tcp_if_include = eno1np0``
    - The network interface that MPI's TCP traffic uses.
    - The setting does, on any node with several interfaces; the interface name is this cluster's
      own. List yours with :command:`ip link`.
  * - ``btl_base_warn_component_unused = 0``
    - Silences a warning about unused network components.
    - Yes, but it is cosmetic. Ubuntu and Debian ship it.

These settings were in place for every measurement and were not compared with Open MPI's defaults,
so no gain in speed is claimed for them. Before copying them, check what your own installation
sets with :command:`ompi_info --all --parsable | grep -E ":param:(pml|btl|btl_tcp_if_include):"`.

To use a setting, put it in the system-wide file, in :file:`~/.openmpi/mca-params.conf`, or on the
command line (``mpirun --mca btl_tcp_if_include <interface> ...``). The file must be readable on
every node: each MPI process reads it itself, and :command:`mpirun` does not pass its contents on.

**Check that a run succeeded from its output.** A run has succeeded only if its output ends with
``Simulation completed successfully``. Neither the job's state nor the exit status is enough: the
script above always ends with ``exit 0``, and :command:`mpirun` returned 0 for a job that lost a
node.

**On Ethernet, check flow control.** If runs become far slower beyond a few nodes (with Open MPI
over TCP, many retransmissions and an MPI test program as slow as FEPX), check that flow control
(IEEE 802.3x) is enabled both on the switch and on the nodes' network cards (:command:`ethtool -a
<interface>`). On the cluster below, with flow control disabled on the switch, FEPX was more than
an order of magnitude slower on 4 nodes, and enabling it was the whole cure.

**Writing to a shared disk during a run slows MPI when both use the same network.** Print only
what you need if a run spends long in its output (the closing summary's *Output*).

**When measuring times,** run one job at a time, since jobs side by side share the network.

.. _performance:

Parallel Performance
--------------------

This section records what to expect from FEPX on more processes, measured on reference cases.

**Measurements on a 4-core laptop** (Intel Core i7-1068NG7, gfortran 15.2, Open MPI 5.0.10, one
process per core). The case: a regularized Voronoi polycrystal (Neper ``-reg 1``) meshed with
10-node tetrahedra (``-order 2 -rcl 0.72``), an FCC material, uniaxial grip loading to a strain
of 0.005 in 2 increments (``initialization isotropic_viscoplastic``, ``midnode_linearization
on``, otherwise the default solver options, e.g. ``cg_preconditioner block_jacobi``), and the
mesh partitioned by Neper for each number of processes.

.. list-table:: Strong scaling: the same mesh on more processes
   :header-rows: 1

   * - Mesh
     - Processes
     - Wall time [s]
     - Speed-up
     - Efficiency
   * - 100 grains, 25,863 elements
     - 1 / 2 / 4
     - 273 / 155 / 92
     - 1 / 1.76 / 2.96
     - 1 / 0.88 / 0.74
   * - 400 grains, 102,014 elements
     - 1 / 2 / 4
     - 1215 / 699 / 447
     - 1 / 1.74 / 2.72
     - 1 / 0.87 / 0.68

.. list-table:: Weak scaling: about 26,000 elements per process
   :header-rows: 1

   * - Mesh
     - Processes
     - Wall time [s]
     - CG iterations per global iteration
     - Efficiency
   * - 100 grains, 25,863 elements
     - 1
     - 273
     - 696
     - 1
   * - 200 grains, 51,812 elements
     - 2
     - 310
     - 883
     - 0.88
   * - 400 grains, 102,014 elements
     - 4
     - 447
     - 1,090
     - 0.61

What these show:

- **The processes do not wait for each other.** In a profiled run (a small mesh on 4
  processes), each process spent 4 to 5 % of its time in communication, and the processes
  did the same amount of work to within 1 %.
- **The element and material work** (the point solves and element matrices; the closing
  summary's *Element and material work*) scales best, slowed only by the processor's clock,
  which falls when all its cores are busy.
- **The conjugate-gradient solves** are limited by memory bandwidth, which the cores of a
  processor share; they scale less well, and take a larger part of the run on larger meshes.
- **On a larger mesh the solver needs more CG iterations per solve** (the third column of the
  weak-scaling table), so the time per process grows even at a constant mesh size per process.
- **One process per core.** Running 8 processes on the 4 cores (hyperthreads) was slower than
  4 processes.

The closing summary of every run (see :ref:`terminal_output`) gives the split of the wall time
between the element and material work and the CG solves, which shows which of the two limits a
given run.

.. _performance_cluster:

Measurements on a Cluster
~~~~~~~~~~~~~~~~~~~~~~~~~

The machine: nodes with one AMD EPYC 7702P (64 cores, 8 memory channels) and 250 GB of memory, on
10 Gb/s Ethernet (TCP), gfortran 15.2, Open MPI 5.0.10, 64 processes per node, one job at a time.
The case: 2000 grains (Neper ``-morpho gg -reg 1``), 10-node tetrahedra, an FCC material with
saturation hardening, uniaxial grip loading to a strain of 0.1 in 100 increments
(``initialization isotropic_viscoplastic``, ``midnode_linearization on``, otherwise the default
solver options), the mesh partitioned by Neper for each number of processes, printing only the
stress. Repeated runs differed by 0.3 to 1.8 %.

.. list-table:: Strong scaling: the same mesh on more nodes
   :header-rows: 1

   * - Mesh
     - Nodes
     - Elements per process
     - Wall time [s]
     - Speed-up
     - Efficiency
   * - 504,469 elements
     - 1
     - 7,882
     - 13,998
     - 1
     - 1
   * -
     - 2
     - 3,941
     - 7,198
     - 1.94
     - 0.97
   * -
     - 4
     - 1,971
     - 3,816
     - 3.67
     - 0.92
   * -
     - 8
     - 985
     - 2,258
     - 6.20
     - 0.78
   * -
     - 12
     - 657
     - 1,335
     - 10.48
     - 0.87
   * -
     - 16
     - 493
     - 849
     - 16.49
     - 1.03
   * - 1,389,557 elements
     - 2
     - 10,856
     - 25,831
     - 1 (from 2 nodes)
     - 1
   * -
     - 4
     - 5,428
     - 13,317
     - 1.94
     - 0.97
   * -
     - 8
     - 2,714
     - 7,586
     - 3.40
     - 0.85
   * -
     - 16
     - 1,357
     - 4,916
     - 5.25
     - 0.66

.. list-table:: Weak scaling: refining the mesh of the same 2000 grains with the nodes
   :header-rows: 1

   * - Mesh
     - Nodes
     - Elements per process
     - Wall time [s]
     - CG iterations
     - Time per CG iteration [ms]
   * - 229,353 elements
     - 1
     - 3,584
     - 5,925
     - 402,652
     - 13.2
   * - 397,470 elements
     - 2
     - 3,105
     - 5,269
     - 406,240
     - 11.8
   * - 723,352 elements
     - 4
     - 2,825
     - 5,853
     - 488,407
     - 11.0
   * - 1,389,557 elements
     - 8
     - 2,714
     - 7,586
     - 602,820
     - 11.8
   * - 2,730,920 elements
     - 16
     - 2,667
     - 11,097
     - 747,382
     - 14.1

What these show, and how to size a job:

- **More nodes were faster on every mesh, up to 16 nodes.** The results here do not extend beyond
  16 nodes, 2.7 million elements or this network. A faster network should help at 8 nodes and
  more.
- **Estimating a run time.** The time per CG iteration is about 3.7 µs per element per process
  (elements per process × 3.7 µs) on 1 or 2 nodes, and grows with the nodes to 4.4 to 4.8 µs on 8 nodes
  and 5.3 to 5.7 µs on 16 nodes (above about 1,000 elements per process). Multiply by the number of CG
  iterations of the run, which depends on the problem far more than on the number of processes.
  CG takes 88 to 95 % of the time on these runs, and the closing summary gives the iterations of a
  short trial run.
- **Run time on refined meshes grows mainly through the CG iterations.** In the weak-scaling table the
  run takes 1.87 times as long on 16 nodes as on 1, and needs 1.86 times as many CG iterations;
  the time per iteration barely changes. The number of processes itself does not change the
  iterations: 443,165 to 443,456 on the 504,469-element mesh on 1 to 16 nodes.
- **Efficiency is not a simple function of the elements per process.** It is lowest (0.66 to 0.85)
  on 8 to 16 nodes at 1,000 to 3,000 elements per process, where communication costs most compared
  with the work saved. It rises again at fewer than about 700 elements per process, which explains the
  efficiency above 1 of the smaller mesh on 16 nodes: do not read it as FEPX beating ideal scaling.
  On one node the time per CG iteration per element is 3.4 to 3.7 µs at 900 to 7,900 elements
  per process but 0.7 µs at 444 (a smaller mesh, 250 grains), and the same 229,353-element
  mesh takes 3.7 µs on 1 node but 2.7 µs on 8 nodes, network included. The likely cause is that each
  process's data then fits in the processor's cache; this was not measured directly.

Restarting a Simulation
-----------------------

A simulation run with ``print restart`` leaves a restart point at the end of every step, and can
be continued with more steps: write a configuration that lists only the steps that follow, add
``restart on``, and run FEPX again in the same directory, on the same number of processes, as
any run. The restart point, what the configuration may change and the output of the continued
run are described in :ref:`restart`; :ref:`tutorial_restart` is a worked example.

.. _terminal_output:

The Terminal Output
-------------------

FEPX writes its progress to the standard output, in the layout of Neper: every line
starts with ``Info   :``, ``Warning:`` or ``Error  :``, and the items are nested by depth, with the
markers ``-``, ``>`` and ``.``; ``[i]`` marks a file read, ``[o]`` a file written. No line is wider
than 72 columns: a longer message continues on the next lines, aligned under its text. Only
process 0 writes (messages of the other processes, such as a warning about a node of their part,
are tagged ``[Rank <r>]`` and written by process 0). The run of the tutorial
:ref:`tutorial_simple` on 4 processes writes (the iterations of increments 3 to 10
left out):

.. code-block:: none

  ========================    F   E   P   X    ===========================
  Info   : A finite element software package for polycrystal plasticity.
  Info   : Version 3.0.0 (git 465464d2, 2026-10-26)
  Info   : Built with: GCC 15.2.0|Open MPI 5.0.10|HDF5 1.14.6
  Info   : Running on 4 ranks.
  Info   : <https://fepx.info>
  Info   : Copyright (C) 1996-2026, and MIT licensed, by the FEPX authors.
  Info   : ---------------------------------------------------------------
  Info   : Loading simulation...
  Info   :   [i] Parsing file `simulation.cfg'...
  Info   :   [i] Parsed file `simulation.cfg'.
  Info   :   [i] Parsing file `simulation.msh'...
  Info   :   [i] Parsed file `simulation.msh'.
  Info   :   - Material parameters:
  Info   :     > Number of phases: 1
  Info   :     > Phase 1: fcc, crystal_plasticity
  Info   :       . c11:                       2.4500E+05 MPa
  Info   :       . c12:                       1.5500E+05 MPa
  Info   :       . c44:                       6.2500E+04 MPa
  Info   :       . m:                         5.0000E-02
  Info   :       . gammadot_0:                1.0000E+00 1/s
  Info   :       . hardening:                 saturation_isotropic
  Info   :       . tau_0:                     2.1000E+02 MPa
  Info   :       . h_0:                       2.0000E+02 MPa
  Info   :       . tau_s:                     3.3000E+02 MPa
  Info   :       . n:                         1.0000E+00
  Info   :       . saturation_evolution:      off (def.)
  Info   :       . cyclic_hardening:          off (def.)
  Info   :       . precipitate_strengthening: off (def.)
  Info   :   - Loading: boundary_conditions general, 2 steps
  Info   :     > Step 1: t = 0.0000E+00 -> 5.0000E+00 s, 5 increments
  Info   :     > Step 2: t = 5.0000E+00 -> 1.0000E+01 s, 5 increments
  Info   :   - Options other than the defaults: (none)
  Info   :   - Output: coo ori stress strain forces (format ascii)
  Info   :   - Mesh parameters:
  Info   :     > Node number: 4008
  Info   :     > Elt  number: 2453 (10-node tetrahedra)
  Info   :     > Elsets: 20; partitions: 4
  Info   :   - Partitioning: 4 ranks
  Info   :     > Elements per rank: 613 to 614
  Info   :     > Ghost nodes per rank: 150 to 191 (645 in total)
  Info   :     > Neighbors per rank: 3 to 3
  Info   :   [o] Writing results to `simulation.sim'...
  Info   : Initializing simulation...
  Info   :   - Elastic step of 1.1351E-02 s (step 1)
  Info   :     > Iteration 1 (SA): R = 2.5329E+00 (613 CG iterations)
  Info   :     > Iteration 2 (SA): R = 1.5028E-05 (612 CG iterations)
  Info   :     > Converged in 2 iterations (1.034 s)
  Info   : Running step 1...
  Info   :   - Increment 1: t = 1.0000E+00 s, dt = 9.8865E-01 s
  Info   :     > Iteration 1 (SA): R = 7.1922E-04 (616 CG iterations)
  Info   :     > Iteration 2 (NR): R = 6.3444E-07 (634 CG iterations)
  Info   :     > Converged in 2 iterations (1.146 s)
  Info   :   - Increment 2: t = 2.0000E+00 s, dt = 1.0000E+00 s
  Info   :     > Iteration 1 (SA): R = 5.1176E-02 (636 CG iterations)
  Info   :     > Iteration 2 (SA): R = 1.7281E-02 (634 CG iterations)
  Info   :     > Iteration 3 (SA): R = 7.3162E-03 (638 CG iterations)
  Info   :     > Iteration 4 (NR): R = 6.5567E-03 (690 CG iterations)
  Info   :     > Iteration 5 (NR): R = 3.8336E-05 (690 CG iterations)
  Info   :     > Converged in 5 iterations (3.188 s)
  ...
  Info   : Final step terminated. Simulation completed successfully.
  Info   :   [o] Wrote results to `simulation.sim'.
  Info   :   - Increments: 10
  Info   :   - Global iterations: 47 (22 NR)
  Info   :   - CG iterations: 36900
  Info   :   - Warnings issued: 1
  Info   :   - Wall time per phase:
  Info   :     > Input and setup:                      0.041 s
  Info   :     > Initialization:                       1.041 s
  Info   :     > Increments:                          34.550 s
  Info   :       . Element and material work:         14.783 s
  Info   :       . CG:                                19.487 s
  Info   :       . Other:                              0.279 s
  Info   :     > Output:                               0.097 s
  Info   : Elapsed time:                              35.728 s
  ========================================================================

The parts of the output:

- **The banner**: the version, with the git commit FEPX was built from (as ``--version``), the
  compiler and the libraries it was built with, and the number of processes.
- **Loading simulation**: the input files as they are read (with the warnings found in them),
  then a summary of the setup: every material parameter of every phase (a parameter not given
  is shown with its default value and ``(def.)``), the loading (boundary conditions,
  constraints, and the plan of every step: its times, increments and strain), the options given
  with a value other than their default (see :ref:`config_file`), the outputs requested, the
  mesh, and its distribution over the processes.
- **Initializing simulation**: the initialization of a crystal-plasticity run (the elastic step,
  or the isotropic viscoplastic start, see :ref:`solver_parameters`).
- **Running step**: each increment, with its end time ``t`` and length ``dt``, then one line per
  global iteration: its method (``SA``, successive approximation, or ``NR``, Newton-Raphson),
  ``R``, the relative change of the velocities that the iteration made (the increment has
  converged when ``R`` is below ``global_tol``), and the conjugate-gradient iterations of its
  linear solve; then the number of iterations and the wall time of the increment. Other events
  of an increment (a rejected attempt, the boundary-velocity iteration of a load-controlled
  face, the landing on a load target, a restart point written) are written under it.
- **The closing summary**: the increments, the global iterations (of them Newton-Raphson) and the
  conjugate-gradient iterations of the run (after the initialization), the number of warnings,
  and the wall time of each phase of the run (input and setup; initialization; the increments,
  split into the element and material work, the conjugate-gradient solves and the rest; the
  output files), then the total elapsed wall time.

Every iteration is also recorded in :file:`simulation.sim/results/convergence` when asked for
with ``print convergence`` (see :ref:`simulation_results`).

**Errors.** A fatal error writes one ``Error  :`` message, ends the run on every process and
leaves a non-zero exit status. The message goes to the standard output (the run's log) and to
the standard error (a batch job's error file), so that it is not lost when one of them is
redirected; in a terminal it therefore shows twice. When several processes stop, every one
writes its message to the standard error (those of processes other than 0 tagged
``[Rank <r>]``), and the standard output has it once. The MPI launcher then reports the end of
the job in its own words. Input errors are found before the simulation starts (see
:ref:`config_file`), and name the file and its line:

.. code-block:: none

  Error  : cfg: simulation.cfg:36: unknown keyword "bogus_keyword"

The results are written to :file:`simulation.sim` (see :ref:`simulation_results`). An existing
:file:`simulation.sim` is not overwritten: FEPX stops with an error, so remove or rename it
first (or continue it from a restart point, see :ref:`restart`).
