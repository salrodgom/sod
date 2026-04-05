set datafile separator comma
set terminal pngcairo size 1800,1200 enhanced font ',10'
set output '/home/salvador/coding/sod/sod/examples/im-18_mc/n04/analysis_plots.png'
set key outside
set grid
set tics out
set border lw 1.2
set multiplot layout 3,2 title 'OUTSOD / ENERGIES analysis: /home/salvador/coding/sod/sod/examples/im-18_mc/n04'

set title 'Pair effect vs distance'
set xlabel 'Ge-Ge distance (A)'
set ylabel 'Delta E vs weighted mean (eV)'
set cblabel 'Weighted count'
set palette rgb 33,13,10
plot '/home/salvador/coding/sod/sod/examples/im-18_mc/n04/pair_energy_report.csv' using 9:8:4 with points pt 7 ps 1.3 lc palette title 'pairs'

set title 'Radial pair curves'
set xlabel 'Ge-Ge distance (A)'
set ylabel 'Pair fraction'
unset cblabel
set ytics nomirror
unset y2tics
plot '/home/salvador/coding/sod/sod/examples/im-18_mc/n04/radial_curve_report.csv' using 3:5 with lines lw 2 title 'Random', \
     '/home/salvador/coding/sod/sod/examples/im-18_mc/n04/radial_curve_report.csv' using 3:7 with lines lw 2 title 'All configs', \
     '/home/salvador/coding/sod/sod/examples/im-18_mc/n04/radial_curve_report.csv' using 3:9 with lines lw 2 title 'Low-energy'

set title 'Radial enrichment'
set xlabel 'Ge-Ge distance (A)'
set ylabel 'Observed / random'
plot '/home/salvador/coding/sod/sod/examples/im-18_mc/n04/radial_curve_report.csv' using 3:10 with lines lw 2 title 'All / random', \
     '/home/salvador/coding/sod/sod/examples/im-18_mc/n04/radial_curve_report.csv' using 3:11 with lines lw 2 title 'Low-energy / random', \
     1.0 with lines dt 2 lw 1 title 'Ratio = 1'

set title 'Radial-bin trends'
set xlabel 'Pair distance (A)'
set ylabel 'Energy trend (eV)'
unset cblabel
set y2label 'Weighted pair fraction'
set y2tics
set ytics nomirror
plot '/home/salvador/coding/sod/sod/examples/im-18_mc/n04/radial_energy_report.csv' using 3:9 with linespoints lw 2 pt 7 title 'Slope per pair', \
     '/home/salvador/coding/sod/sod/examples/im-18_mc/n04/radial_energy_report.csv' using 3:8 with linespoints lw 2 pt 5 title 'Presence delta', \
     '/home/salvador/coding/sod/sod/examples/im-18_mc/n04/radial_energy_report.csv' using 3:6 axes x1y2 with impulses lw 2 title 'Pair fraction'
unset y2tics
set ytics auto
unset y2label

set title 'Cluster slopes'
set xlabel 'Cluster cutoff distance (A)'
set ylabel 'Slope per cluster (eV)'
plot '/home/salvador/coding/sod/sod/examples/im-18_mc/n04/cluster_energy_report.csv' using 1:($2==2 ? $6 : 1/0) with linespoints lw 2 pt 7 title 'size 2', \
     '/home/salvador/coding/sod/sod/examples/im-18_mc/n04/cluster_energy_report.csv' using 1:($2==3 ? $6 : 1/0) with linespoints lw 2 pt 7 title 'size 3', \
     '/home/salvador/coding/sod/sod/examples/im-18_mc/n04/cluster_energy_report.csv' using 1:($2==4 ? $6 : 1/0) with linespoints lw 2 pt 7 title 'size 4'

set title 'Cluster presence delta'
set xlabel 'Cluster cutoff distance (A)'
set ylabel 'Presence delta (eV)'
plot '/home/salvador/coding/sod/sod/examples/im-18_mc/n04/cluster_energy_report.csv' using 1:($2==2 ? $5 : 1/0) with linespoints lw 2 pt 7 title 'size 2', \
     '/home/salvador/coding/sod/sod/examples/im-18_mc/n04/cluster_energy_report.csv' using 1:($2==3 ? $5 : 1/0) with linespoints lw 2 pt 7 title 'size 3', \
     '/home/salvador/coding/sod/sod/examples/im-18_mc/n04/cluster_energy_report.csv' using 1:($2==4 ? $5 : 1/0) with linespoints lw 2 pt 7 title 'size 4'

unset multiplot
