#!/bin/bash

# Post-process the `Uniaxial Control' example (`1_uniaxial')
# - This bash script requires a full installation of Neper v4.2.0
# - This bash script requires gnuplot

# Merge *.core* files and write *.step* files
neper -S .

# Visualize the 33 component of elemental stress
neper -V ../1_uniaxial.sim \
    -simstep 2 \
    -datanodecoo coo \
    -dataelt1drad 0.004 \
    -dataelt3dedgerad 0.0015 \
    -datanodecoofact 10 \
    -dataelt3dcol stress33 \
    -dataeltscaletitle "Stress 33 (MPa)" \
    -dataeltscale 0:850 \
    -showelt1d all \
    -cameraangle 13.5 \
    -imagesize 800:800 \
    -print 1_s33_deform

# Visualize the 33 component of elemental strain
neper -V ../1_uniaxial.sim \
    -simstep 2 \
    -datanodecoo coo \
    -dataelt1drad 0.004 \
    -dataelt3dedgerad 0.0015 \
    -datanodecoofact 10 \
    -dataelt3dcol strain-el33 \
    -dataeltscaletitle "Strain 33 (-)" \
    -dataeltscale 0.000:0.001:0.002:0.003:0.004 \
    -showelt1d all \
    -cameraangle 13.5 \
    -imagesize 800:800 \
    -print 1_e33_deform

# Plot the stress-strain curve with gnuplot
gnuplot plot_stress_strain.gp

# This script produces the files:
# - 1_s33_deform.png
# - 1_s33_deform-scale3d.png
# - 1_e33_deform.png
# - 1_e33_deform-scale3d.png
# - 1_stress_strain.png

exit 0
