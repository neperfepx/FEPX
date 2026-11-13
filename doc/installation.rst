.. _installation:

Installation
============

FEPX can be installed from a conda package (Linux and macOS) or built from its source with
CMake. Both give the executable :file:`fepx`, its manual page (``man fepx``) and the tutorial
cases of :ref:`tutorials`, ready to run.

.. _installation_conda:

Conda Package
-------------

FEPX is available as a conda package from `conda-forge <https://conda-forge.org>`_, for Linux and
macOS, with the HDF5 output, in two builds: with Open MPI and with MPICH. No compilation and no
administrator rights are needed.

If conda is not available on your system (command ``conda``), install it first, using
`Miniforge <https://github.com/conda-forge/miniforge>`_:

.. code-block:: console

  $ curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
  $ bash Miniforge3-$(uname)-$(uname -m).sh

Accept the license and the default installation directory (:file:`$HOME/miniforge3`), and answer
``yes`` to the last question, about initializing conda, so that ``conda`` is available in new
terminals. Then, open a new terminal; the installer file can be removed.

Then, install FEPX in a new environment (named ``fepx`` here) and activate it:

.. code-block:: console

  $ conda create -n fepx -c conda-forge fepx
  $ conda activate fepx

This installs the Open MPI build, with Open MPI; for the MPICH build, use ``"fepx=*=mpich*"``
instead of ``fepx``. FEPX is then run as described in :ref:`running`, with the ``mpirun`` of the
environment, and the tutorial cases are in :file:`$CONDA_PREFIX/share/doc/fepx/examples/` (see
:ref:`installation_examples`). Update FEPX with ``conda update -n fepx -c conda-forge fepx``.

FEPX and `Neper <https://neper.info>`_ can also be installed together, in the same environment
(named ``neperfepx`` here), to generate meshes with Neper, simulate them with FEPX and
post-process the results with Neper:

.. code-block:: console

  $ conda create -n neperfepx -c conda-forge neper fepx
  $ conda activate neperfepx

The MPI of a conda environment is meant for a workstation: on a computer cluster, build FEPX from
source with the cluster's MPI (see :ref:`installation_source`), so that it uses the cluster's
network.

.. _installation_source:

Building from Source
--------------------

Requirements
~~~~~~~~~~~~

FEPX is written in Fortran 2018 and runs on Unix-like systems. It needs:

- a Fortran 2018 compiler: GNU gfortran (tested) or Intel ifx (generally supported, though not tested);
- an MPI implementation providing the ``mpi_f08`` Fortran module (Open MPI or MPICH);
- CMake 3.20 or newer;
- optionally, for the HDF5 output (``output_format hdf5``, see :ref:`simulation_results`): the
  serial HDF5 library with its Fortran interface, and for its tests the HDF5 tools; a
  thread-safe HDF5 built with CMake (e.g. conda-forge's) also needs a C compiler (e.g. gcc,
  which comes with gfortran), to be found;
- for the test suite: Python 3 with NumPy.

On Debian and Ubuntu:

.. code-block:: console

  $ sudo apt install gfortran libopenmpi-dev openmpi-bin cmake python3-numpy
  $ sudo apt install libhdf5-dev hdf5-tools    # for the HDF5 output

On macOS, with Homebrew:

.. code-block:: console

  $ brew install gcc open-mpi cmake numpy
  $ brew install hdf5-mpi                      # for the HDF5 output

Homebrew's ``hdf5-mpi`` is a parallel HDF5 built with Homebrew's Open MPI, which FEPX then uses
serially (the configuration warns about it); it works because FEPX is built with the same Open
MPI. On arm64 processors (Apple silicon, and aarch64 Linux), FEPX is compiled without fused
multiply-add contraction (``-ffp-contract=off``), so that its results are those of an x86-64
build.

The Source
~~~~~~~~~~

Each release provides a source tarball, :file:`fepx-<version>.tar.gz`, and a separate test-data
tarball, :file:`fepx-<version>-testdata.tar.gz`, which holds the reference results of the
regression tests (see :ref:`installation_testing`). Unpack the source, and, to run the
regression tests, the test data over it (both unpack into :file:`fepx-<version>/`):

.. code-block:: console

  $ tar -xzf fepx-<version>.tar.gz
  $ tar -xzf fepx-<version>-testdata.tar.gz      # optional: the regression tests
  $ cd fepx-<version>

Build and Install
~~~~~~~~~~~~~~~~~

From the top directory of the source:

.. code-block:: console

  $ cmake -S src -B build -DCMAKE_BUILD_TYPE=Release
  $ cmake --build build -j
  $ cmake --build build --target userinstall

This builds :file:`build/fepx` and installs it for the user (no administrator rights needed),
under :file:`$HOME/.local`, in the layout of the conda package: :file:`bin/fepx`,
:file:`share/man/man1/fepx.1.gz` and :file:`share/doc/fepx/` (the last command is also ``make
userinstall`` in :file:`build`, as for Neper). Add :file:`$HOME/.local/bin` to the ``PATH``
if it is not there already (a warning says so), e.g. in :file:`~/.bashrc`:

.. code-block:: console

  export PATH="$HOME/.local/bin:$PATH"

To install for all users instead, under :file:`/usr/local` (or the prefix given by
``-DCMAKE_INSTALL_PREFIX``), run ``sudo cmake --install build`` (``sudo make install`` in
:file:`build`). The build directory can also be used without installing: run
:file:`build/fepx`.

CMake Options
~~~~~~~~~~~~~

Options are given to the first command, ``cmake -S src -B build``, as ``-D<option>=<value>``:

- ``-DCMAKE_BUILD_TYPE=<type>``: how the code is compiled, which sets its speed and the checks
  made while it runs:

  - ``Release`` (the default): optimized, for production runs;
  - ``Debug``: no optimization, with run-time checks (array bounds, floating-point exceptions,
    uninitialized values), for use with a debugger; 4 to 6 times slower than ``Release``;
  - ``Checked``: the run-time checks of ``Debug`` in optimized code, to look for errors in
    longer runs; about 1.5 times slower than ``Release``.

- ``-DCMAKE_INSTALL_PREFIX=<directory>``: where ``cmake --install`` installs FEPX (default
  :file:`/usr/local`; the target ``userinstall`` always installs under :file:`$HOME/.local`).

- ``-DFEPX_ENABLE_HDF5=<value>``: the optional HDF5 output (``output_format hdf5``):

  - ``AUTO`` (the default): built when the serial HDF5 library with its Fortran interface is
    found (and, with the test suite, the HDF5 tools ``h5ls``, ``h5dump``, ``h5diff``; without
    them a warning says that the HDF5 output is not built);
  - ``ON``: required (the configuration fails without them);
  - ``OFF``: FEPX is built without HDF5, writes the text files only, and ``output_format hdf5``
    is then an error that says so.

  The configuration's last line says which: ``HDF5 output: ON`` or ``OFF``. A serial HDF5 is
  used (only one process writes the file), so any MPI works with it.

- ``-DMPI_Fortran_COMPILER=<MPI Fortran wrapper>``: the MPI to build with, when several are
  installed (default: the one CMake finds first, usually ``mpif90``). For example, MPICH on
  Debian and Ubuntu: ``-DMPI_Fortran_COMPILER=mpif90.mpich
  -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec.mpich``. FEPX must then be run with that MPI's
  ``mpirun``.

- ``-DFEPX_INSTALL_RPATH=<value>``: whether the installed executable keeps the paths of the MPI
  and HDF5 libraries it was built with:

  - ``ON`` (the default): it runs without ``LD_LIBRARY_PATH`` when they are outside the system
    directories (e.g. from an environment module or :file:`/opt`);
  - ``OFF``: for packages (e.g. the conda packages, whose run paths conda sets).

- ``-DFEPX_WERROR=<value>``: ``ON`` turns compiler warnings into errors (default ``OFF``).

- ``-DBUILD_TESTING=<value>``: ``OFF`` skips the test suite and its Python requirement (default
  ``ON``).

The conda packages are made from the recipe of :file:`conda/recipe` with `rattler-build
<https://rattler-build.prefix.dev>`_: ``make conda`` in the build directory builds them from the
source tree (one per MPI implementation) and runs their tests, into :file:`conda/<platform>/`.

.. _installation_testing:

Testing
~~~~~~~

FEPX comes with unit tests, verification tests (analytic solutions), regression tests (against
earlier FEPX results) and MPI tests on 1, 2 and 4 ranks. Run them from the top directory with:

.. code-block:: console

  $ ctest --test-dir build --output-on-failure

Labels select subsets, e.g. ``ctest --test-dir build -L mpi`` (MPI tests only), ``-L
regression`` or ``-L docs`` (the documentation's checks). Three labels divide the suite by run time:

- ``quick``: the unit tests and a selection of the other tests that reaches every source file
  and every configuration keyword, a few minutes in all (``ctest --test-dir build -L quick``);
- ``full``: every test but the tutorial cases, the default check of a change;
- ``tutorial``: the tutorial cases (about 30 minutes), which ``full`` leaves out.

The regression tests need the reference results of the test-data tarball; without them they are
not defined, and the configuration says so.

.. _installation_examples:

The Installed Tutorial Cases
----------------------------

:file:`share/doc/fepx/examples/` of the installation (:file:`$CONDA_PREFIX/share/doc/fepx/examples/`
for the conda package) holds one directory per tutorial case of :ref:`tutorials`, with its configuration file, its mesh as
:file:`simulation.msh` and its other input files. To run one, copy its directory and run FEPX
in the copy:

.. code-block:: console

  $ cp -r $CONDA_PREFIX/share/doc/fepx/examples/simple_simulation .
  $ cd simple_simulation
  $ mpirun -np 4 fepx

The case :file:`restart_cycle2` continues the simulation of :file:`restart_cycle1` (see
:ref:`restart`), and :file:`postprocessing/` holds the script of the post-processing tutorial.
