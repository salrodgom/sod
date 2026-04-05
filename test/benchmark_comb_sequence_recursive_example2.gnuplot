set terminal pngcairo size 1200,800 enhanced font 'Helvetica,12'
set output '/home/salvador/coding/sod/sod/test/benchmark_comb_sequence_recursive_example2.png'
set datafile separator ','
set key left top
set xlabel 'N substitutions'
set ylabel 'Wall time per level in sequential sweep (s)'
set title 'example2 sequential sweep from clean dir: recursion-enabled comb comparison'
set grid ytics xtics
set logscale y
set xtics 1
plot \
  '/home/salvador/coding/sod/sod/test/benchmark_comb_sequence_recursive_example2.csv' using 1:2 with linespoints lw 2 pt 5 ps 1.2 title 'classic 0.62 combsod', \
  '/home/salvador/coding/sod/sod/test/benchmark_comb_sequence_recursive_example2.csv' using 1:3 with linespoints lw 2 pt 9 ps 1.2 title 'sod_ensemble comb'
