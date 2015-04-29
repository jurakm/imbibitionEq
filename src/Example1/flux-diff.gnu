set xlabel "time"
set ylabel "Nonwetting source"
set key  center top
set grid
set xrange [0:1]
 set terminal postscript eps color solid lw 3
# set output "diff_Q.eps"
 set output "diff_Q_bv.eps"

plot "clin-flux.txt" u 1:2   t "constant linearized",\
     "vlin-flux.txt" u 1:3   t "variable linearized",\
     "nlin-flux.txt" u 1:2 w l t "nonlinear",\
     "nlin-flux.txt" u 1:($4 -0.35)*0.05 w l t "boundary value"


pause -1
