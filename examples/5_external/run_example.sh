#!/bin/bash

# Run the `Running a Simulation with External Definition Files' example
# (`5_external')
# - Executes FEPX with OpenMPI's `mpirun' on 4 processors by default.

# Run FEPX in parallel on 4 processors
mpirun -np 4 fepx

exit 0

# This script produces the files:
# - post.report
# - post.coo.core*
# - post.ori.core*
# - post.defrate-eq.core*
# - post.force.*
# - post.work.core*
