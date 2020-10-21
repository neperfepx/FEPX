#!/bin/bash

# Generate the tessellation and mesh files for use in example `4_restart'
# - This bash script requires a configured installation of Neper 4.0.0
#   in order to be properly executed.

# First, generate a centroidal tessellation for the domain with `Neper -T':
neper -T -n 100 -reg 1 -rsel 1.25 -mloop 4 \
    -morpho "diameq:1,1-sphericity:lognormal(0.145,0.03)" -morphooptistop val=5e-3 \
    -oricrysym "cubic" \
    -o simulation

# Then, generate a coarse finite element mesh for the domain with `Neper -M':
neper -M simulation.tess -order 2 -rcl 1.25 -part 2

# Copy the mesh file into the first simulation directory for execution:
cp simulation.msh cycle1/.

exit 0

# This script produces output files:
# - simuation.tess
# - simulation.msh
