.. _introduction:

Introduction
============

Description
-----------

FEPX is a finite element software package for polycrystal plasticity. It is well-suited to model the global and local mechanical behavior of large polycrystalline solids as aggregates of grains as well as associated microstructural evolution through a highly scalable parallel framework. Each grain is discretized by finite elements whose local behavior corresponds accordingly to the local behavior of a sub-volume of a crystal. These behaviors include:

- Nonlinear kinematics capable of resolving large (or finite) strains and large rotations,

- Anisotropic elasticity based on crystal symmetry,

- Anisotropic plasticity based on rate-dependent slip restricted to dominant slip systems,

- Appropriate state variable evolution for crystal lattice orientation and slip system strengths.

FEPX strives to be a user-friendly, efficient, and robust tool. All of the input data are prescribed non-interactively, using ASCII files. FEPX works hand in hand with `Neper <https://neper.info>`_, which can be used for both generating the input polycrystal mesh and post-processing the simulation results. Currently, FEPX is designed to provide solutions on rectangular prismatic domains (or right cuboids).

Resources and Support
---------------------

Several complementary resources describing FEPX are available:

- The FEPX reference manual, the document you are reading, provides a detailed overview of all of FEPX's capabilities. Specific sections are dedicated to simulation input and output, running simulations, and various example simulations.

- The FEPX theory manual, written by Paul Dawson and Donald Boyce, provides in depth details on the underlying mechanical theory and finite element methods utilized in FEPX. It is available at `<https://arxiv.org/abs/1504.03296>`_ (please note that the descriptions of simulation input and output provided in the FEPX theory manual are no longer up-to-date and the user is instead recommended to utilize the descriptions provided in the FEPX reference manual).

- The FEPX website, `<https://fepx.info>`_, gives a general introduction to FEPX with illustrative examples.

- The `FEPX GitHub repository <https://github.com/acmelab-ua/FEPX>`_ is where the latest version is available and where user interactions take place:

  - To get and keep up-to-date with the latest version, clone the repository using:

    .. code-block:: console

      $ git clone https://github.com/acmelab-ua/FEPX.git

    which gives access to the latest stable development release on the default, :code:`main` branch. To update your local repository, run :command:`git pull` from within the repository.

  - To report bugs, use the `issue tracker <https://github.com/acmelab-ua/FEPX/issues>`_. When reporting bugs, please provide a minimal working example and the FEPX terminal output.

  - To ask questions, share comments or request new features, use `discussions <https://github.com/acmelab-ua/FEPX/discussions>`_.

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

- Install FEPX on your system (as root):

  .. code-block:: console

    $ make install

This procedure uses the default configuration options and should work out-of-the-box if you have a Fortran compiler, OpenMPI, and CMake installed. Testing is performed on GFortran version 6 and greater, and OpenMPI version 2 and greater (other Fortran compilers and MPI distributions may also work, though they are not explicitly supported or tested by ACME Lab). A minimum version of CMake version 3.0 is required to utilize the build system.

Testing FEPX
------------

FEPX comes packaged with tests and reference outputs. To run the tests, execute the following from your build folder:

.. code-block:: console

  $ make test

or (equivalently):

.. code-block:: console

  $ ctest

This runs the tests in :code:`Normal` mode, for which the produced output files are compared to reference ones.

The (packaged) reference output files are generated on Ubuntu version 20.04, using compiler GFortran version 9.3.0} and OpenMPI version 4.0.3 (note: CMake will switch to the MPI Fortran compiler to build FEPX, which will be built against GFortran version 9.3.0). It is expected that different versions may result in minor (insignificant) changes to the output, though this will generally result in failed tests. If this happens, you may switch to the :code:`Minimal` mode as described in the following.

The testing mode is controlled by variable :code:`BUILD_TESTING_MODE`, which may be changed using:

.. code-block:: console

  $ ccmake ..

for an interactive command-line tool, or:

.. code-block:: console

  $ cmake-gui ..

for an interactive graphical tool, or directly at the command line, using Cmake's :data:`-D` option:

.. code-block:: console

  $ cmake -DBUILD_TESTING_MODE={Normal,Minimal,Writing} ..

- The (default) :code:`Normal` mode checks if the program completes without error and if the produced output is the same as a set of reference output.

- The :code:`Minimal` mode only checks if the program completes without error. This mode may be useful when installing on a machine which has program or library versions different from the ones with which the reference output was generated.

- The :code:`Writing` mode overwrites the reference outputs with the generated output.  This mode may be useful when installing on a machine which has program or library versions different from the ones with which the reference output was generated and the user needs a reference output before making changes to the source code.

Getting Started
---------------

To run a serial simulation on a local computer, the :program:`fepx` binary must be run in a terminal:

.. code-block:: console

  $ fepx

or, for parallel simulations:

.. code-block:: console

  $ mpirun -np <N> fepx

where :data:`<N>` refers to the number of MPI processes (typically equal to or less than the number of cores on the local machine). The :program:`fepx` binary should always be run from within a simulation directory that contains the necessary simulation input files (:ref:`Simulation Input <simulation_input>`).

To perform simulations across multiple computational nodes on an HPC cluster, a submission script that conforms to the specific job scheduling program is necessary. Examples of generic scripts for common job scheduling programs are detailed in :ref:`Running a Simulation <running_a_simulation>`.

During a simulation run, FEPX returns real-time messages in the terminal and, upon successful completion, prints requested output data in ASCII files.

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

The development of FEPX began in the late 1990s and was lead by Paul Dawson, and involved many members of the Deformation Process Laboratory (DPLab) at Cornell University, until early 2020.  An extended development history contributed by Paul Dawson, the lead investigator of the DPLab, can be found in :ref:`Development History <development_history>`.

Ongoing development has since been lead by Matthew Kasemer, and involved other members of the Advanced Computational Materials Engineering Laboratory (ACME Lab) at The University of Alabama.
