#!/bin/bash

# Run the `Multiaxial Control with Constant Load Rate' example (`3_triaxCLR')
# - Executes FEPX with OpenMPI's `mpirun' on 4 processors by default.

# Run FEPX in parallel on 4 processors
mpirun -np 4 fepx

exit 0

# This script produces the files:
# - post.report
# - post.coo.core*
# - post.crss.core*
# - post.strain_eq.core*
# - post.force.*
