#!/bin/bash

# Post-process the `Multiaxial Control with Constant Strain Rate' example (`2_triaxCSR') in Neper
# - This bash script will generate the associated figures in the documentation
#   for this example. The post-processed simulation directory `2_triaxCSR.sim'
#   is placed in the `examples/' directory.
#
# - This bash script requires a configured installation of Neper 4.0.0
#   in order to be properly executed.

# First, the completed simulation directory is parsed through `Neper -S'
neper -S .

# Next, use `Neper -V' to visualize the elemental equivalent plastic strain
neper --rcfile none -V ../2_triaxCSR.sim -simstep 3 -datanodecoo coo \
        -dataelt1drad 0.004 -dataelt3dedgerad 0.0015 -datanodecoofact 10 \
        -dataelt3dcol strain-pl-eq -showelt1d all \
        -dataeltscale 0.000:0.010 -dataeltscaletitle "Equivalent Plastic Strain" \
        -cameraangle 13.5 -imagesize 800:800 \
        -print 2_eqplstrain_preunload

# Then, use `Neper -V' to visualize the elemental plastic work
neper --rcfile none -V ../2_triaxCSR.sim -simstep 3 -datanodecoo coo \
        -dataelt1drad 0.004 -dataelt3dedgerad 0.0015 -datanodecoofact 10 \
        -dataelt3dcol work-pl -dataeltscaletitle "Plastic Work (mJ)" \
        -dataeltscale 0.00:6.50 -showelt1d all \
        -cameraangle 13.5 -imagesize 800:800 -print 2_plwork_preunload

# Finally, create the strain-time curve with gnuplot
sed -i '3i 0 0 0 0 0 1 0 1' post.force.x1
sed -i '3i 0 0 0 0 0 1 0 1' post.force.y1
sed -i '3i 0 0 0 0 0 1 0 1' post.force.z1
gnuplot gen_norm_straintime_plot.gp

# This script produces output files:
# - 2_eqplstrain_preunload.png
# - 2_eqplstrain_preunload-scale3d.png
# - 2_plwork_preunload.png
# - 2_plwork_preunload-scale3d.png
# - 2_normalstraintime.png

exit 0
