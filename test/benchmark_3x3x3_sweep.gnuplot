set terminal pngcairo size 900,600 enhanced font 'Sans,13'
set output '/home/salvador/coding/sod/sod/test/benchmark_3x3x3_sweep.png'

set title "SOD comb enumeration: Modern vs Legacy (quartz 3×3×3, 81 sites)" font 'Sans,15'
set xlabel "Substitution level N"
set ylabel "Wall time (s)"

set datafile separator ','
set key top left box opaque
set grid ytics lc rgb '#dddddd'
set style data linespoints

set logscale y
set format y "%.3g"
set yrange [0.01:*]
set xrange [-0.3:6.3]
set xtics 1

# Add a secondary y-axis label for speedup
set y2label "Speedup (legacy / modern)"
set y2tics
set y2range [1:2.5]
set format y2 "%.1f"

plot '/home/salvador/coding/sod/sod/test/benchmark_3x3x3_sweep.csv' \
     every ::1 using 1:4 with linespoints lw 2 pt 7 ps 1.3 lc rgb '#2166ac' title 'Modern SOD', \
     '' every ::1 using 1:5 with linespoints lw 2 pt 5 ps 1.3 lc rgb '#b2182b' title 'Legacy combsod', \
     '' every ::1 using 1:($5/$4) axes x1y2 with linespoints lw 1.5 pt 9 ps 1.2 lc rgb '#4daf4a' dt 2 title 'Speedup (right axis)'
