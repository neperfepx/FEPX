.. _introduction:

Introduction
============

Description
-----------

FEPX is a finite element software package for polycrystal plasticity. It is well-suited to model the global and local mechanical behavior of large polycrystalline solids as aggregates of grains as well as associated microstructural evolution through a highly scalable parallel framework. Each grain is discretized by finite elements which have a crystal behavior. This includes:

- Nonlinear kinematics capable of resolving large (or finite) strains and large rotations,

- Anisotropic elasticity based on crystal symmetry,

- Anisotropic plasticity based on rate-dependent slip restricted to dominant slip systems,

- Appropriate state variable evolution for crystal lattice orientation and slip system strengths.

FEPX strives to be a user-friendly, efficient, and robust tool. All of the input data are prescribed non-interactively, using ASCII files. FEPX works hand in hand with `Neper <https://neper.info>`_, which can be used for both generating the input polycrystal mesh and post-processing the simulation results. As a finite element code, FEPX allows for complex loadings on general domains, but it also offers capabilities specific to RVE-type calculations (uniform loading of cubic polycrystals).

Resources and Support
---------------------

Several complementary resources describing FEPX are available:

- The FEPX reference manual, the document you are reading, provides a detailed overview of all of FEPX's capabilities. Specific sections are dedicated to simulation input and output, running simulations, and various example simulations.

- The FEPX website, `<https://fepx.info>`_, gives a general introduction to FEPX with illustrative examples.

- The `FEPX GitHub repository <https://github.com/neperfepx/fepx>`_ is where the latest version is available and where user interactions take place:

  - To get and keep up-to-date with the latest version, clone the repository using:

    .. code-block:: console

      $ git clone https://github.com/neperfepx/fepx.git

    which gives access to the latest stable development release on the default, :code:`main` branch. To update your local repository, run :command:`git pull` from within the repository.

  - To report bugs, use the `issue tracker <https://github.com/neperfepx/fepx/issues>`_. When reporting bugs, please provide a minimal working example and the FEPX terminal output.

  - To ask questions, share comments or request new features, use `discussions <https://github.com/neperfepx/fepx/discussions>`_.

- The two FEPX reference papers:

  - `R Quey and M Kasemer, The Neper/FEPX project: free / open-source polycrystal generation, deformation simulation, and post-processing, IOP Conference Series Materials Science and Engineering, vol. 1249, pp. 012021, 2022 <https://iopscience.iop.org/article/10.1088/1757-899X/1249/1/012021/meta>`_.

  - `P.R. Dawson and D.E. Boyce, FEPX -- Finite Element Polycrystals: Theory, Finite Element Formulation, Numerical Implementation and Illustrative Examples, arXiv, 2015 <https://arxiv.org/abs/1504.03296>`_ describes the underlying mechanical theory and finite element methods (please note that the provided descriptions of simulation input and output are no longer up-to-date).

Resources for Neper can be accessed from `<https://neper.info>`_.

Installing FEPX
---------------

FEPX is written in Fortran, and it can run on any Unix-like system (including MacOS). Parallelization of the code is achieved via OpenMPI. Compilation is performed via CMake:

- Create a :file:`build` directory, for instance as a subdirectory of FEPX's :file:`src` directory:

  .. code-block:: console

    $ cd src
    $ mkdir build

- Run CMake from within the :file:`build` directory, pointing to FEPX's :file:`src` directory:

  .. code-block:: console

    $ cd build
    $ cmake ..

- Build FEPX:

  .. code-block:: console

    $ make

  This generates the :file:`fepx` binary file.  To make it available system-wide, run (as root):

  .. code-block:: console

    $ make install

This procedure uses the default configuration options and should work out-of-the-box if you have a Fortran compiler, OpenMPI, and CMake installed. Testing is performed on GFortran version 6 and greater, and OpenMPI version 2 and greater (other Fortran compilers and MPI distributions may also work, though they are not explicitly supported or tested). A minimum version of CMake version 3.0 is required to utilize the build system.

Fine Configuration
~~~~~~~~~~~~~~~~~~

Program configuration options can be specified, before installation, using CMake variables.   This can be done using a terminal tool:

.. code-block:: console

  $ ccmake ..

or an interactive tool:

.. code-block:: console

  $ cmake-gui ..

or directly at the command line, using Cmake's :data:`-D` option:

.. code-block:: console

  $ cmake -D<VARIABLE1>=<VALUE1> -D<VARIABLE2>=<VALUE2> ..

The program configuration variable concerns the precision of the output:

- :code:`OUTPUT_PRECISION`, default :data:`1e-12` (12 decimal digits).

  Values are printed to this precision, and smaller values (in absolute value) are printed as :data:`0`.

Testing FEPX
------------

FEPX comes packaged with tests and reference outputs. To run the tests, execute the following from your build folder:

.. code-block:: console

  $ make test

or (equivalently):

.. code-block:: console

  $ ctest

This runs the tests in :code:`Normal` mode, for which the produced output files are compared to reference ones.

The (packaged) reference output files are generated on Ubuntu version 22.04, using compiler GFortran version 11.4.0 and OpenMPI version 4.1.2 (note: CMake will switch to the MPI Fortran compiler to build FEPX, which will be built against GFortran version 11.4.0). It is expected that different versions may result in minor (insignificant) changes to the output, though this will generally result in failed tests. If this happens, you may switch to the :code:`Minimal` mode as described in the following. The testing mode is controlled by variable :code:`BUILD_TESTING_MODE`, which may be changed as described in :ref:`fine_configuration`.

- The (default) :code:`Normal` mode checks if the program completes without error and if the produced output is the same as a set of reference output.

- The :code:`Minimal` mode only checks if the program completes without error. This mode may be useful when installing on a machine which has program or library versions different from the ones with which the reference output was generated.

- The :code:`Writing` mode overwrites the reference outputs with the generated output.  This mode may be useful when installing on a machine which has program or library versions different from the ones with which the reference output was generated and the user needs a reference output before making changes to the source code.

In the standard case where the testing mode is :code:`Normal`, it is further possible to choose how the produced output is compared to the reference output using the variable :code:`BUILD_TESTING_DIFF`, which can take the value of :code:`Soft` to allow for small differences (the default, recommended) or :code:`Hard` otherwise.

.. note:: As the program is tested on a different environment from the one on which the reference files were generated (specified above), it may happen that some tests fail.  This problem can generally be ignored.

Getting Started
---------------

Running a simulation requires two files: a :ref:`mesh_file`, which describes the material, and a configuration file, which defines the simulation itselt, as described in :ref:`config_file` and related sections.  With these files at disposal, to run a serial simulation on a local computer, the :program:`fepx` binary must be run in a terminal (from the directory where the input files are located):

.. code-block:: console

  $ fepx

For parallel simulations, FEPX must be run as:

.. code-block:: console

  $ mpirun -np <N> fepx

where :data:`<N>` refers to the number of MPI processes (typically equal to or less than the number of cores on the local machine).

To perform simulations across multiple computational nodes on an HPC cluster, a submission script that conforms to the specific job scheduling program is necessary (see examples in :ref:`hpc`).

During a simulation run, FEPX returns real-time messages in the terminal and, upon successful completion, prints requested output data in a :ref:`simulation_directory` named :file:`simulation.sim`.  The :file:`.sim` format is a human-friendly database format which is shared with (and readable by) :ref:`Neper <https://neper.info>`.

Conventions Used in This Manual
-------------------------------

- A command entered at the terminal is shown like this:

  .. code-block:: console

    $ command

  The first character on the line is the terminal prompt, and should not be typed. The dollar sign, :data:`$`, is used as the standard prompt in this manual, although some systems may use a different character.

- A program (or command) option is printed like :data:`this`;

- An option argument is printed like :data:`<this>`;

- The name of a variable, or a meta-syntactic variable, is printed like :data:`<this>`;

- Literal examples are printed like

  .. code-block:: console

    this

- File names are printed like :file:`this`.

Additionally, hereinafter a :data:`core` will explicitly refer to a processor (or CPU) of a computer. This terminology is also consistent with file name formatting for parallel simulation output by FEPX.

Development History
-------------------

The development of FEPX began in the late 1990s and was lead by Paul Dawson, and involved many members of the Deformation Process Laboratory (DPLab) at Cornell University, USA, until early 2020.  An extended development history contributed by Paul Dawson, the lead investigator of the DPLab, can be found in :ref:`development_history`.  Ongoing development has since been lead by Matthew Kasemer (Advanced Computational Materials Engineering Laboratory (ACME Lab), The University of Alabama, USA), and Romain Quey (CNRS, Mines Saint-Étienne, France), and involved other members of their respective groups.
