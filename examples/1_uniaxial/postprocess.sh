#!/bin/bash

# Post-process the `Uniaxial Control' example (`1_uniaxial') in Neper
# - This bash script will generate the associated figures in the documentation
#   for this example. The post-processed simulation directory `1_uniaxial.sim'
#   is placed in the `examples/' directory.
#
# - This bash script requires a configured installation of Neper 4.0.0
#   in order to be properly executed.

# First, the completed simulation directory is parsed through `Neper -S'
neper -S .

# Next, use `Neper -V' to visualize the 33 component of elemental stress
neper --rcfile none -V ../1_uniaxial.sim -simstep 4 -datanodecoo coo \
        -dataelt1drad 0.004 -dataelt3dedgerad 0.0015 -datanodecoofact 10 \
        -dataelt3dcol stress33 -dataeltscaletitle "Stress 33 (MPa)" \
        -dataeltscale 0:1000 -showelt1d all \
        -cameraangle 13.5 -imagesize 800:800 -print 1_s33_deform

# Then, use `Neper -V' to visualize the 33 component of elemental strain
neper --rcfile none -V ../1_uniaxial.sim -simstep 4 -datanodecoo coo \
        -dataelt1drad 0.004 -dataelt3dedgerad 0.0015 -datanodecoofact 10 \
        -dataelt3dcol strain-el33 -dataeltscaletitle "Strain 33 (-)" \
        -dataeltscale 0.000:0.001:0.002:0.003:0.004 -showelt1d all \
        -cameraangle 13.5 -imagesize 800:800 -print 1_e33_deform

# Finally, create the stress-strain curve with gnuplot
sed -i '3i 0 0 0 0 0 1 0' post.force.z1
gnuplot gen_stressstrain_plot.gp

# This script produces output files in `examples/1_uniaxial':
# - 1_s33_deform.png
# - 1_s33_deform-scale3d.png
# - 1_e33_deform.png
# - 1_e33_deform-scale3d.png
# - 1_stressstrain.png

exit 0
