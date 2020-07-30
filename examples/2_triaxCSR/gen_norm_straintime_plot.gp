# FEPX Post-Process Gnuplot Script
# This script does NOT pad the initial zeros
#    instead the postprocess.sh script adds these via `sed'.

set term pngcairo enhanced font 'Latin Modern Roman,18'
set output '2_normalstraintime.png'

set key out vert right center box width 1
set grid xtics ytics dt 3 lw 1 lc -1

set xlabel "Time (s)"
set ylabel "Normal Strains (-)"
set size ratio 0.75
plot "post.force.x1" using ($7):($8 - (1)) w lp lc 1 lw 2 pt 5  ps 1 title "E11", \
     "post.force.y1" using ($7):($8 - (1)) w lp lc 2 lw 2 pt 7  ps 1 title "E22", \
     "post.force.z1" using ($7):($8 - (1)) w lp lc 3 lw 2 pt 13 ps 1 title "E33"
