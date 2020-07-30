#!/bin/bash

# Visualize the polycrystal and finite element mesh in Neper
# - This bash script will generate the associated images for Figure 5.1.
# - This bash script requires a configured installation of Neper 4.0.0
#   in order to be properly executed.

# First, the morphology of the polycrystal
neper --rcfile none -V example_3.tess -dataedgerad 0.004 -datacellcol ori \
      -imagesize 800:800 -cameraangle 13.5 \
      -print 0_morpho

# Then, the finite element mesh overlaid on the morphology
neper --rcfile none -V example_3.tess,simulation.msh \
      -dataelt1drad 0.004 -dataelt3dedgerad 0.0015 \
      -showelt1d all -imagesize 800:800 \
      -cameraangle 13.5 -dataeltcol ori \
      -print 0_femesh

# This script produces output files in `mesh/':
# - 0_morpho.png
# - 0_femesh.png

exit 0
