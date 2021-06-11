#!/bin/bash

# Generate the tessellation and mesh files for the `Uniaxial Control' example
# (`1_uniaxial')
# - This bash script requires a full installation of Neper v4.2.0

# Generate a tessellation describing a microstructure geometry
neper -T \
    -n 50 \
    -reg 1 \
    -morpho gg \
    -o simulation

# Generate a finite element mesh
neper -M simulation.tess \
    -order 2 \
    -part 4

exit 0

# This script produces the files:
# - simuation.tess
# - simulation.msh
