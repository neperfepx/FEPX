#!/bin/bash
#
# Post-process the `Running a Simulation with External Definition Files'
#   example (`5_external') in Neper
#
# - This bash script will generate the associated figures in the documentation
#   for this example. The post-processed simulation directory `2_triaxCSR.sim'
#   is placed in the `examples/' directory.
#
# - This bash script requires a configured installation of Neper 4.0.0
#   in order to be properly executed.
#
# First, the completed simulation directory is parsed through `Neper -S'
neper -S .

#
# Next, use `Neper -V' to visualize the elemental equivalent deformation rate
neper --rcfile none -V ../5_external.sim -simstep 1 -datanodecoo coo \
        -dataelt1drad 0.004 -dataelt3dedgerad 0.0015 -datanodecoofact 50 \
        -dataelt3dcol defrate-eq -dataeltscaletitle "Equivalent Deformation Rate" \
        -dataeltscale 0.000:0.030 -showelt1d all \
        -imagesize 800:800 -cameraangle 13.5 -print 5_defrateeq
#
# Use `Neper -V' to visualize the elemental total work
neper --rcfile none -V ../5_external.sim -simstep 2 -datanodecoo coo \
        -dataelt1drad 0.004 -dataelt3dedgerad 0.0015 -datanodecoofact 50 \
        -dataelt3dcol work -dataeltscaletitle "Work (mJ)" \
        -dataeltscale 0.00:6.50 -showelt1d all \
        -imagesize 800:800 -cameraangle 13.5 -print 5_work
#
# Optionally, visualize the per-element orientations from the mesh
# - Note: The per-element orientations need to be embedded into the mesh file
#   by running `generate_mesh.sh'.
#./generate_mesh.sh
# neper --rcfile none -V simulation.msh -dataelt3dedgerad 0 -dataelt1drad 0.004 \
#     -showelt1d all \
#     -imagesize 800:800 -cameraangle 13.5 \
#     -dataeltcol ori \
#     -print 5_eltorientations

# Finally, create the stress-time curve with gnuplot
sed -i '3i 0 0 0 0 0 1 0' post.force.z1
gnuplot gen_stresstime_plot.gp

#
# This script produces output files:
# - 5_defrateeq.png
# - 5_defrateeq-scale3d.png
# - 5_work.png
# - 5_work-scale3d.png
# - 5_eltorientations.png (OPTIONAL)
# - 5_stresstime.png

exit 0
