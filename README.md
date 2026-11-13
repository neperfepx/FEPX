# FEPX - Finite Element Polycrystal Plasticity

## Overview

FEPX is a parallel finite element program for elasto-viscoplastic simulations of the large
deformations of polycrystalline aggregates, written in modern Fortran (2018) with MPI. Its
primary features are:

- Nonlinear kinematics, for large strains and large rotations
- Anisotropic elasticity based on crystal symmetry (fcc, bcc, hcp and bct crystals)
- Anisotropic plasticity based on rate-dependent slip on the slip systems of the crystal
- Hardening models: saturation (isotropic or anisotropic, with optional cyclic hardening and
  precipitate strengthening) and Kocks-Mecking dislocation density (with optional irradiation
  softening)
- State variable evolution for the crystal lattice orientations and the slip system strengths
- Loading by boundary-condition presets, general velocity conditions, multi-point constraints
  or periodic conditions, with strain, time, load and stress targets
- 10-node tetrahedral and 20-node hexahedral elements
- Results independent of the number of MPI processes
- Output as a simulation directory read by Neper, or in HDF5 with an XDMF file for ParaView

FEPX works with [Neper](https://neper.info), which generates its input meshes and post-processes
and visualizes its results.

## Installation

FEPX is available as a conda package from [conda-forge](https://conda-forge.org), for Linux and
macOS:

```
conda create -n fepx -c conda-forge fepx
```

FEPX can also be built from source, which requires a Fortran 2018 compiler (e.g. GNU gfortran),
an MPI implementation providing the `mpi_f08` module (Open MPI or MPICH) and CMake (3.20 or
later); the tests require Python 3 with NumPy, and the optional HDF5 output the serial HDF5
library with its Fortran interface:

```
cmake -S src -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build -L quick
cmake --install build --prefix <dir>
```

The installation instructions, including the build options, are in the
[documentation](https://fepx.info/doc/index.html) (`doc/installation.rst`).

## Documentation

The documentation, with tutorials, is available at <https://fepx.info/doc/index.html> and in
`doc/`; `man fepx` summarizes the program's usage.

## License

MIT, see [`COPYING`](COPYING).
