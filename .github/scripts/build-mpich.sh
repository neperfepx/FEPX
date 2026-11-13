#!/usr/bin/env bash
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Build MPICH from source for the CI jobs (ROADMAP step 4.3).
#
# Ubuntu 24.04's MPICH package links libmpich against PMIx while its mpiexec (Hydra) provides
# only PMI-1, so every process of "mpiexec -n N" starts as a singleton (1 rank); the MPI tests
# then ran N independent copies. This build uses MPICH's own Hydra launcher and PMI, and the
# ch3:nemesis device (shared memory and TCP, no libfabric or UCX), which is all single-node
# tests need. The install is cached by the workflow (key: version and this script).
#
# Usage: build-mpich.sh VERSION PREFIX
#   VERSION  MPICH release (e.g. 4.2.3), downloaded from mpich.org
#   PREFIX   install directory (bin/mpifort, bin/mpiexec, ...)

set -euo pipefail

version="$1"
prefix="$2"
work="$(mktemp -d)"

# Download and unpack the release.
curl --fail --location --retry 5 --retry-delay 10 --silent --show-error \
  --output "${work}/mpich.tar.gz" \
  "https://www.mpich.org/static/downloads/${version}/mpich-${version}.tar.gz"
tar -xzf "${work}/mpich.tar.gz" -C "${work}"

# Configure, build and install. gfortran >= 10 needs -fallow-argument-mismatch for MPICH's
# legacy Fortran bindings (the mpi_f08 module that FEPX uses is unaffected).
cd "${work}/mpich-${version}"
./configure --prefix="${prefix}" \
  --with-device=ch3:nemesis --with-pm=hydra \
  --enable-fortran=all --disable-cxx --disable-doc \
  CC=gcc FC=gfortran F77=gfortran \
  FFLAGS=-fallow-argument-mismatch FCFLAGS=-fallow-argument-mismatch \
  > "${work}/configure.log" 2>&1 || { tail -50 "${work}/configure.log"; exit 1; }
make -j "$(nproc)" > "${work}/make.log" 2>&1 || { tail -50 "${work}/make.log"; exit 1; }
make install > "${work}/install.log" 2>&1 || { tail -50 "${work}/install.log"; exit 1; }

# Report the installed version.
"${prefix}/bin/mpichversion" | head -3
rm -rf "${work}"
