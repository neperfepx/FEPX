#!/bin/bash

# Post-process the `Multiaxial Control with Constant Load Rate' example
# (`3_triaxCLR')
# - This bash script requires a full installation of Neper v4.2.0
# - This bash script requires gnuplot

# Merge *.core* files and write *.step* files
neper -S .

# Visualize the elemental critical resolved shear stress on the first slip system
neper -V ../3_triaxCLR.sim \
    -simstep 4 \
    -datanodecoo coo \
    -dataelt1drad 0.004 \
    -dataelt3dedgerad 0.0015 \
    -dataelt3dcol crss1 \
    -dataeltscaletitle "Critical Resolved Shear Stress (MPa)" \
    -showelt1d all \
    -dataeltscale 200.0:202.0 \
    -cameraangle 13.5 \
    -imagesize 800:800 \
    -print 3_crss_deform

# Visualize the grain/phase distribution from tessellation
neper -V simulation.tess \
    -datacellcol group \
    -cameraangle 13.5 \
    -imagesize 800:800 \
    -print 3_phase

# Plot the stress-time and stress-strain curves with gnuplot
gnuplot plot_stress_strain.gp
gnuplot plot_stress_time.gp

# This script produces the files:
# - 3_crss_deform.png
# - 3_crss_deform-scale3d.png
# - 3_phase.png
# - 3_stress_time.png
# - 3_stress_strain_33.png

exit 0
