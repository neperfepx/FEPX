#!/bin/bash

NEPER="neper --rcfile none"

# a more simple microstructure
#N=20
#RCL=1
#DOMAIN="cube(0.25,0.25,1)"
#MORPHO="gg"

# exact microstructure on website
N="from_morpho"
RCL=0.6
DOMAIN="cube(0.3,0.15,1.25)"
MORPHO="diameq:lognormal(0.035,0.01)+0.2*lognormal(0.075,0.02),1-sphericity:lognormal(0.145,0.03)"

# tessellation

$NEPER -T -n $N \
          -reg 1 -rsel $RCL \
          -domain $DOMAIN \
          -morpho $MORPHO \
          -morphooptistop "val=1e-2,xreps=0.001" \
          -group "mode!=2?1:2" \
          -o specimen

# grain and phase visualizations

$NEPER -V specimen.tess \
          -datacellcol id \
          -datacelltrs 0.25 \
          -dataedgerad 0.002 \
          -cameracoo 0.125+vx:0.125+vy:0.5+vz \
          -cameraangle 4.5 \
          -imagesize 1000:3000 \
          -imageformat pov,png \
          -print specimen-tess

$NEPER -V specimen.tess \
          -datacellcol group \
          -datacelltrs 0.25 \
          -dataedgerad 0.002 \
          -cameracoo 0.125+vx:0.125+vy:0.5+vz \
          -cameraangle 4.5 \
          -imagesize 1000:3000 \
          -imageformat pov,png \
          -print specimen-phase

# meshing

$NEPER -M specimen.tess -rcl $RCL -order 2 -part 4:16

# grain and phase mesh visualizations

$NEPER -V specimen.tess,specimen.msh \
          -dataelsetcol id \
 	  -dataelt1drad 0.0015 \
          -dataelt3dedgerad 0.0003 \
          -dataelt3dedgecol 16:16:16 \
          -showelt1d all \
          -cameracoo 0.125+vx:0.125+vy:0.5+vz \
          -cameraangle 4.5 \
          -imagesize 1000:3000 \
          -imageformat pov,png \
          -print specimen-mesh
         
$NEPER -V specimen.tess,specimen.msh \
          -dataelsetcol group \
          -dataelt1drad 0.0015 \
          -dataelt3dedgerad 0.0003 \
          -dataelt3dedgecol 16:16:16 \
          -showelt1d all \
          -cameracoo 0.125+vx:0.125+vy:0.5+vz \
          -cameraangle 4.5 \
          -imagesize 1000:3000 \
          -imageformat pov,png \
          -print specimen-mesh-phase

# simulation

cp specimen.tess simulation/simulation.tess
cp specimen.msh simulation/simulation.msh
cd simulation
mpirun -np 8 fepx

# post-processing

cd ..
$NEPER -S simulation

# field variables visualizations

$NEPER -V simulation.sim \
          -cameracoo 0.125+vx:0.125+vy:0.5+vz \
          -cameralookat 0.125:0.125:0.5 \
          -cameraangle 4.5 \
          -dataelt1drad 0.0015 \
          -dataelt3dedgerad 0.0000 \
          -dataelt3dedgecol 16:16:16 \
          -showelt1d all \
          -imagesize 1000:3000 \
          -dataeltscale 0.00:0.1 \
          -dataeltcolscheme viridis \
          -dataeltscaletitle "Equivalent strain" \
          -loop STEP 0 1 7 \
            -simstep STEP \
            -datanodecoo from_sim \
            -dataeltcol strain-eq \
            -print strain-eq-stepSTEP \
          -endloop

$NEPER -V simulation.sim \
          -cameracoo 0.125+vx:0.125+vy:0.5+vz \
          -cameralookat 0.125:0.125:0.5 \
          -cameraangle 4.5 \
          -dataelt1drad 0.0015 \
          -dataelt3dedgerad 0.0000 \
          -dataelt3dedgecol 16:16:16 \
          -showelt1d all \
          -imagesize 1000:3000 \
          -dataeltscale 200:900 \
          -dataeltcolscheme viridis \
          -dataeltscaletitle "Stress33" \
          -loop STEP 0 1 7 \
            -simstep STEP \
            -datanodecoo from_sim \
            -dataeltcol stress33 \
            -print stress33-stepSTEP \
          -endloop

exit 0
