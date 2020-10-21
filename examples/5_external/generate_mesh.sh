#!/bin/bash

# Generate the tessellation and mesh files for use in example `5_external'
# - This bash script requires a configured installation of Neper 4.0.0
#   in order to be properly executed.

# First, copy over the `orifile' from an external directory
cp ../mesh/orifile .

# Next, generate a dual-phase centroidal tessellation for the domain with `Neper -T':
neper -T -n 100 -reg 1 -rsel 1.25 -mloop 4 \
    -morpho "diameq:1,1-sphericity:lognormal(0.145,0.03)" -morphooptistop val=5e-3 \
    -ori "file(orifile,des=rodrigues)" -oridistrib "normal(5)" \
    -oridescriptor "rodrigues:active" \
    -o simulation

# Then, generate a coarse finite element mesh for the domain with `Neper -M':
neper -M simulation.tess -order 2 -rcl 1.25 \
    -format "msh,ori" -part 2 \
    -o simulation

# Remove the copied `orifile'
rm orifile

exit 0

# This script produces output files:
# - simuation.tess
# - simulation.msh
# - simulation.ori
