# FEPX Post-Process Gnuplot Script
# This script does NOT pad initial zeros

set term pngcairo enhanced font 'Latin Modern Roman,18'
set output '5_stresstime.png'

unset key
set grid xtics ytics dt 3 lw 1 lc -1

set xlabel "Time (s)"
set ylabel "Stress (MPa)"
plot "post.force.z1" using ($7):($5 / $6) w lp lc -1 lw 2 pt 13 ps 1

