set term pngcairo enhanced font 'Latin Modern Roman,15'
set output 'stress_strain.png'

unset key
set grid xtics ytics lw 1

set xlabel "Strain [%]"
set ylabel "Stress [MPa]"
set xrange [0:1]
set yrange [0:400]
set size ratio 0.75
plot "post.force.z1" using ($7 * 0.1):($5) w lp lc 3 lw 1.5 ps 1.5
