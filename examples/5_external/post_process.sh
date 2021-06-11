#!/bin/bash
#
# Post-process the `Running a Simulation with External Definition Files'
# example (`5_external')
# - This bash script requires a full installation of Neper v4.2.0
# - This bash script requires gnuplot

# Merge *.core* files and write *.step* files
neper -S .

# Visualize the elemental equivalent deformation rate
neper -V ../5_external.sim \
    -simstep 1 \
    -datanodecoo coo \
    -dataelt1drad 0.004 \
    -dataelt3dedgerad 0.0015 \
    -datanodecoofact 50 \
    -dataelt3dcol defrate-eq \
    -dataeltscaletitle "Equivalent Deformation Rate" \
    -dataeltscale 0.000:0.030 \
    -showelt1d all \
    -imagesize 800:800 \
    -cameraangle 13.5 \
    -print 5_defrate_eq

# Visualize the elemental total work
neper -V ../5_external.sim \
    -simstep 2 \
    -datanodecoo coo \
    -dataelt1drad 0.004 \
    -dataelt3dedgerad 0.0015 \
    -datanodecoofact 50 \
    -dataelt3dcol work \
    -dataeltscaletitle "Work (mJ)" \
    -dataeltscale 0.00:2.50 \
    -showelt1d all \
    -imagesize 800:800 \
    -cameraangle 13.5 \
    -print 5_work

# Visualize the per-element orientations from the mesh
neper -V simulation.msh \
    -dataelt3dedgerad 0 \
    -dataelt1drad 0.004 \
    -showelt1d all \
    -imagesize 800:800 \
    -cameraangle 13.5 \
    -dataeltcol ori \
    -print 5_elt_orientations

# Plot the stress-time curve with gnuplot
gnuplot plot_stress_time.gp

# This script produces the files:
# - 5_defrate_eq.png
# - 5_defrate_eq-scale3d.png
# - 5_work.png
# - 5_work-scale3d.png
# - 5_elt_orientations.png (OPTIONAL)
# - 5_stress_time.png

exit 0
