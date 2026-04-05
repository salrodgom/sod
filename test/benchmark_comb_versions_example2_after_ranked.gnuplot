set terminal pngcairo size 1200,800 enhanced font 'Helvetica,12'
set output '/home/salvador/coding/sod/sod/test/benchmark_comb_versions_example2_after_ranked.png'
set datafile separator ','
set key left top
set xlabel 'N substitutions'
set ylabel 'Wall time (s)'
set title 'example2 benchmark after ranked visited-core port'
set grid ytics xtics
set logscale y
set xtics 1
plot \
  '/home/salvador/coding/sod/sod/test/benchmark_comb_versions_example2_after_ranked.csv' using 1:2 with linespoints lw 2 pt 7 ps 1.2 title 'legacy local combsod', \
  '/home/salvador/coding/sod/sod/test/benchmark_comb_versions_example2_after_ranked.csv' using 1:3 with linespoints lw 2 pt 5 ps 1.2 title 'classic 0.62 combsod', \
  '/home/salvador/coding/sod/sod/test/benchmark_comb_versions_example2_after_ranked.csv' using 1:4 with linespoints lw 2 pt 9 ps 1.2 title 'sod_ensemble comb'
