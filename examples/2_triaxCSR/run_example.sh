#!/bin/bash

# Run the `Multiaxial Control with Constant Strain Rate' example (`2_triaxCSR')
# - Executes FEPX with OpenMPI's `mpirun' on 4 processors by default.

# Run FEPX in parallel on 4 processors
mpirun -np 4 fepx

exit 0

# This script produces the files:
# - post.report
# - post.coo.core*
# - post.strain-pl-eq.core*
# - post.work-pl.core*
# - post.force.*
