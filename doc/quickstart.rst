.. _quickstart:

Quick Start
===========

A first simulation in five minutes: a polycrystal of 20 grains pulled to 1 % strain, the first
tutorial (:ref:`tutorial_simple`), whose input files are installed with FEPX.

1. Install FEPX (:ref:`installation`), and Neper (`<https://neper.info>`_) to visualize the
   results.

2. Copy the case and go to its directory. It holds the mesh, :file:`simulation.msh` (made by
   Neper: 10-node tetrahedra, one element set per grain), and the configuration,
   :file:`simulation.cfg` (installed in :file:`share/doc/fepx/examples/` of the installation:
   :file:`$CONDA_PREFIX/share/doc/fepx/examples/` for the conda package, see
   :ref:`installation_examples`):

   .. code-block:: console

     $ cp -r $CONDA_PREFIX/share/doc/fepx/examples/simple_simulation .
     $ cd simple_simulation

3. Read the configuration: an fcc material (elastic constants, slip kinetics, saturation
   hardening), velocity conditions on the faces, two steps to 5 s and 10 s, and the variables
   to write:

   .. literalinclude:: ../tutorials/simple_simulation/simulation.cfg

4. Run FEPX on 4 processes (the mesh has 8 partitions: FEPX says so and distributes it evenly):

   .. code-block:: console

     $ mpirun -np 4 fepx

   .. literalinclude:: tutorials/generated/simple_simulation/terminal.txt
      :language: none

5. The results are in the simulation directory, :file:`simulation.sim`:

   .. literalinclude:: tutorials/generated/simple_simulation/tree.txt
      :language: none

   Neper reads it: for example, the volume averages of the stress and the strain at the end of
   each step are written by

   .. literalinclude:: tutorials/generated/simple_simulation/commands.txt
      :language: console
      :lines: 1

   and the stress field of the last step, on the deformed mesh, is plotted with ``neper -V``
   (the commands are in :ref:`tutorial_simple`):

   .. figure:: tutorials/generated/simple_simulation/stress33_step2.png
      :align: center
      :width: 50%

      Stress :math:`\sigma_{33}` at 1 % strain.

Next: :ref:`config_file` describes the configuration file, :ref:`running` the runs, and
:ref:`tutorials` gives more examples.
