set xlabel "time"
set ylabel "Nonwetting source"
set key  left   top
set grid
set xrange [0:1]
set title "delta = 0.05"
set terminal postscript eps color solid lw 3
set output "diff_Q_bv_Ex2-0.05.eps"

plot "del_0.05clin-flux.txt" u 1:2 w l  t "constant linearized",\
     "del_0.05vlin-flux.txt" u 1:3 w l  t "variable linearized",\
     "del_0.05nlin-flux.txt" u 1:2 w l  t "nonlinear",\
     "del_0.05nlin-flux.txt" u 1:($4/100)  t "boundary value/100"

pause -1
