#!/bin/bash

# Generate the tessellation and mesh files for the test cases
# - This bash script requires a full installation of Neper v4.2.0

# Generate a tessellation describing a microstructure geometry
neper -T \
    -n 2 \
    -reg 1 \
    -morphooptiini 'coo:file(seed_points.txt)' \
    -o simulation

# Generate a finite element mesh
neper -M simulation.tess \
    -order 2 \
    -part 2 \
    -rcl 3.5

# Generate visualizations
neper -V simulation.tess \
    -print test_tess
neper -V simulation.msh \
    -print test_mesh

exit 0

# This script produces the files:
# - simuation.tess
# - simulation.msh
