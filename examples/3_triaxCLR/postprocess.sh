#!/bin/bash

# Post-process the `Multiaxial Control with Constant Load Rate' example (`3_triaxCLR') in Neper
# - This bash script will generate the associated figures in the documentation
#   for this example. The post-processed simulation directory `3_triaxCLR.sim'
#   is placed in the `examples/' directory.
#
# - This bash script requires a configured installation of Neper 4.0.0
#   in order to be properly executed.

# First, the completed simulation directory is parsed through `Neper -S'
neper -S .

# Next, use `Neper -V' to visualize the elemental critical resolved shear stress on the first system
neper --rcfile none -V ../3_triaxCLR.sim -simstep 5 -datanodecoo coo \
        -dataelt1drad 0.004 -dataelt3dedgerad 0.0015 \
        -dataelt3dcol crss1 -dataeltscaletitle "Critical Resolved Shear Stress (MPa)" \
        -showelt1d all -dataeltscale 237.0:241.50 \
        -cameraangle 13.5 -imagesize 800:800 \
        -print 3_crss_deform

# Then, use `Neper -V' to visualize the grain/phase distribution from tessellation
cp ../mesh/example_3.tess .

neper --rcfile none -V example_3.tess -datacellcol group \
        -cameraangle 13.5 -imagesize 800:800 -print 3_phase

rm example_3.tess

# Finally, create the stress-time and stress-strain curves with gnuplot
sed -i '3i 0 0 0 0 0 1 0 1' post.force.x1
sed -i '3i 0 0 0 0 0 1 0 1' post.force.y1
sed -i '3i 0 0 0 0 0 1 0 1' post.force.z1
gnuplot gen_norm_stressstrain_plot.gp
gnuplot gen_norm_stresstime_plot.gp

# This script produces output files:
# - 3_crss_deform.png
# - 3_crss_deform-scale3d.png
# - 3_phase.png
# - 3_stresstime.png
# - 3_stressstrain_33.png

exit 0
