#!/bin/bash

# Generate the tessellation and mesh files for the `Running a Simulation with
# External Definition Files' example (`5_external')
# - This bash script requires a full installation of Neper v4.2.0

# Generate a tessellation describing a microstructure geometry
# Load a set of orientations from file, "orifile"
neper -T \
    -n 50 \
    -reg 1 \
    -group "id<=25?1:2" \
    -ori "file(orifile,des=rodrigues)" \
    -oridistrib "normal(5)" \
    -oridescriptor "rodrigues:active" \
    -morpho gg \
    -o simulation

# Generate a finite element mesh
neper -M simulation.tess \
    -order 2 \
    -format "msh,ori" \
    -part 4

exit 0

# This script produces the files:
# - simuation.tess
# - simulation.msh
# - simulation.ori
