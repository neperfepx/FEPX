# FEPX Post-Process Gnuplot Script

set term pngcairo enhanced font 'Latin Modern Roman,18'
set output '3_stress_time.png'

set key out vert right center box width 1
set grid xtics ytics dt 3 lw 1 lc -1

set xlabel "Time (s)"
set ylabel "Normal Stresses (MPa)"

plot "post.force.x1" using ($7):($3 / $6) w lp lc 1 lw 2 pt 5  ps 0.5 title "S11", \
     "post.force.y1" using ($7):($4 / $6) w lp lc 2 lw 2 pt 7  ps 0.5 title "S22", \
     "post.force.z1" using ($7):($5 / $6) w lp lc 3 lw 2 pt 13 ps 0.5 title "S33"
