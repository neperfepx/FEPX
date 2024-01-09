# FEPX Post-Process Gnuplot Script

set term pngcairo enhanced font 'Latin Modern Roman,18'
set output '3_stress_strain_33.png'

unset key
set grid xtics ytics dt 3 lw 1 lc -1

set xlabel "Strain (-)"
set ylabel "Stress (MPa)"
set size ratio 0.75
plot "post.force.z1" using ($8 - 1):($5 / $6) w lp lc -1 lw 2 pt 13 ps 1 title ""
