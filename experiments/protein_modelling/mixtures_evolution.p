set terminal post eps color enhanced
set autoscale	# scale axes automatically
set xtic auto	# set xtics automatically
set ytic auto	# set ytics automatically
set style data linespoints

set xtics font "Times-Roman, 20"
set ytics font "Times-Roman, 20"

set ytics nomirror

set xlabel "Iterations\n"
set ylabel "Number of inferred components\n"
set xlabel font "Times-Roman, 25"
set ylabel font "Times-Roman, 25"
set y2label font "Times-Roman, 25"

set key at 45,15 font ",20" spacing 1.5 

set output "mixtures_evolution.eps"

plot "bvm_sine_mml" using 1:2 title "BVM Sine" lc rgb "blue"
     
