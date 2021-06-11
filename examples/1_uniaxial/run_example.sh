#!/bin/bash

# Run the `Uniaxial Control' example (`1_uniaxial')
# - Executes FEPX with OpenMPI's `mpirun' on 4 processors by default.

# Run FEPX in parallel on 4 processors
mpirun -np 4 fepx

exit 0

# This script produces the files:
# - post.report
# - post.coo.core*
# - post.strain-el.core*
# - post.stress.core*
# - post.force.*
