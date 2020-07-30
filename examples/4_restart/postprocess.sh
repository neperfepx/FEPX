#!/bin/bash

# Post-process the `Multiaxial Control with Constant Strain Rate' example (`4_restart') in Neper
# - This bash script will generate the associated figures in the documentation
#   for this example. The post-processed simulation directory `cycle1.sim' and `cycle2.sim'
#   is placed in the `examples/4_restart/' directory.
#
# - This bash script requires a configured installation of Neper 4.0.0
#   in order to be properly executed.

# First, the completed simulation directory is parsed through `Neper -S'
neper -S cycle1
cp cycle1/post.report cycle2/.
sed -i 's/number_of_steps 4/number_of_steps 8/g' cycle2/post.report
neper -S cycle2

# Next, use `Neper -V' to visualize the elemental accumulated slip shear on (111)[1-10] system:
neper --rcfile none -V cycle2.sim -simstep 8 -datanodecoo coo \
        -dataelt1drad 0.004 -dataelt3dedgerad 0.0015 \
        -dataeltscaletitle "Accumulated Slip Shear" \
        -showelt1d all -dataeltscale -0.030:0.030 \
        -cameraangle 13.5 -imagesize 800:800 \
        -dataelt3dcol slip3 -print 4_slipshear3_cycle2 \
        -dataelt3dcol slip7 -print 4_slipshear7_cycle2

# This script produces output files:
# - 4_slipshear3_cycle2.png
# - 4_slipshear3_cycle2-scale3d.png
# - 4_slipshear7_cycle2.png
# - 4_slipshear7_cycle2-scale3d.png

exit 0
