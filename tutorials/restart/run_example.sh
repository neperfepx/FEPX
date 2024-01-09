#!/bin/bash

# Run the `Restarting a Simulation with Appended Load Steps' example
# (`4_restart')
# - Executes FEPX with OpenMPI's `mpirun' on 4 processors by default.

# Create a directory for the first cycle and copy simulation files
mkdir cycle1
cp simulation_cycle1.config ./cycle1/simulation.config
cp simulation.tess simulation.msh ./cycle1
cd cycle1

# Run FEPX in parallel on 4 processors (first cycle)
mpirun -np 4 fepx

# Create a directory for the second cycle and move simulation files
cd ../
mkdir cycle2
cp simulation_cycle2.config ./cycle2/simulation.config
cp simulation.tess simulation.msh ./cycle2
cp cycle1/rst0.* ./cycle2
cd cycle2

# Run FEPX in parallel on 4 processors (second cycle)
mpirun -np 4 fepx

exit 0

# This script produces the files:
# - post.report
# - post.coo.core*
# - post.defrate_pl_eq.core*
# - post.slip.core*
# - post.sliprate.core*
# - post.restart.control*
# - post.restart.field.*
