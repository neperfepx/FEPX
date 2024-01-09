#!/bin/bash

# Post-process the `Multiaxial Control with Constant Strain Rate' example
# (`2_triaxCSR')
# - This bash script requires a full installation of Neper v4.2.0
# - This bash script requires gnuplot

# Merge *.core* files and write *.step* files
neper -S .

# Visualize the elemental equivalent plastic strain
neper -V ../2_triaxCSR.sim \
    -step 2 \
    -datanodecoo coo \
    -dataelt1drad 0.004 \
    -dataelt3dedgerad 0.0015 \
    -datanodecoofact 10 \
    -dataelt3dcol strain_pl_eq \
    -showelt1d all \
    -dataeltscale 0.000:0.002:0.004:0.006:0.008 \
    -dataeltscaletitle "Equivalent Plastic Strain" \
    -cameraangle 13.5 \
    -imagesize 800:800 \
    -print 2_eqplstrain

# Visualize the elemental plastic work
neper -V ../2_triaxCSR.sim \
    -step 2 \
    -datanodecoo coo \
    -dataelt1drad 0.004 \
    -dataelt3dedgerad 0.0015 \
    -datanodecoofact 10 \
    -dataelt3dcol work_pl \
    -dataeltscaletitle "Plastic Work (mJ)" \
    -dataeltscale 0.00:4.00 \
    -showelt1d all \
    -cameraangle 13.5 \
    -imagesize 800:800 \
    -print 2_plwork

# Plot the strain-time curve with gnuplot
gnuplot plot_strain_time.gp

# This script produces the files:
# - 2_eqplstrain.png
# - 2_eqplstrain-scale3d.png
# - 2_plwork.png
# - 2_plwork-scale3d.png
# - 2_strain_time.png

exit 0
