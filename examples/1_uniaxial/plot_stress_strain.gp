# FEPX Post-Process Gnuplot Script
# This script does NOT pad the initial zeros
#    instead the postprocess.sh script adds these via `sed'.

set term pngcairo enhanced font 'Latin Modern Roman,18'
set output '1_stress_strain.png'

unset key
set grid xtics ytics dt 3 lw 1 lc -1

set xlabel "Strain (-)"
set ylabel "Stress (MPa)"
set xrange [0:0.02]
set yrange [0:550]
set size ratio 0.75
plot "post.force.z1" using ($7 * (1e-2)):($5) w lp lc -1 lw 2 pt 13 ps 1.5
