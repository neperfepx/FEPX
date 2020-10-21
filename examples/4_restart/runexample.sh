#!/bin/bash

# Run the `Restarting a Simulation with Appended Load Steps' example (`4_restart') in FEPX
# - This bash script will copy the necessary mesh file from the external
#   directory `../mesh' by default. The mesh file may also be generated
#   using the provided `generate_mesh.sh' script if Neper is installed.
#
# - Executes FEPX with OpenMPI's `mpirun' on 2 processors by default.

# First, enter the directory of the first cycle
cd cycle1

# Next, copy the mesh file into the directory
cp ../../mesh/simulation.msh .

# Then, run FEPX in parallel
mpirun -np 2 fepx

# Now, enter the directory of the second cycle
cd ../cycle2

# First, copy over simulation data needed to restart the simulation with new files
cp ../cycle1/cycle1.control .
cp ../cycle1/cycle1.field.* .
cp ../cycle1/simulation.msh .

# Then, run FEPX in parallel
mpirun -np 2 fepx

exit 0

# This script produces output files:
# - post.report
# - post.coo.*
# - post.defrate-pl-eq.*
# - post.slip.*
# - post.sliprate.*
