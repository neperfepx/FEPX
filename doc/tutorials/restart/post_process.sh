#!/bin/bash

# Post-process the `Multiaxial Control with Constant Strain Rate' example
# (`4_restart')
# - This bash script requires a full installation of Neper v4.2.0

# Merge *.core* files, write *.step* files, and merge cycles together
neper -S cycle1
neper -S cycle2
neper -S cycle1.sim,cycle2.sim -o cycle1-2.sim

# Visualize the elemental accumulated slip shear on various systems
# during the 2nd step of the 2nd cycle (4th overall step)
neper -V cycle1-2.sim \
    -step 4 \
    -datanodecoo coo \
    -dataelt1drad 0.004 \
    -dataelt3dedgerad 0.0015 \
    -dataeltscaletitle "Accumulated Slip Shear" \
    -showelt1d all \
    -dataeltscale -0.0010:0.0010 \
    -cameraangle 13.5 \
    -imagesize 800:800 \
    -dataelt3dcol slip7 \
    -print 4_slipshear7_step4 \
    -dataelt3dcol slip8 \
    -print 4_slipshear8_step4

# This script produces the files:
# - 4_slipshear7_step4.png
# - 4_slipshear7_step4-scale3d.png
# - 4_slipshear8_step4.png
# - 4_slipshear8_step4-scale3d.png

exit 0
